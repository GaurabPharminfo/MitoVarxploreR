library(shiny)
library(bslib)
library(dplyr)
library(DT)
library(htmltools)
library(tidyr)
library(stringr)

allVariantsUI <- function(id) {
  ns <- NS(id)
  
  layout_sidebar(
    sidebar = sidebar(
      title = "Filters & Tools",
      width = 300,
      
      p(strong("Views")),
      actionButton(ns("view_default"), "Default Table", class = "btn-outline-primary mb-2"),
      actionButton(ns("view_summary"), "Summarise", class = "btn-outline-info mb-2"),
      actionButton(ns("view_breakup"), "Patient Breakup", class = "btn-outline-primary"),
      hr(),
      
      textInput(ns("search_gene"), "Search Gene:", placeholder = "e.g., MT-ND1"),
      textInput(ns("search_variant"), "Search Variant ID(s):", placeholder = "e.g., m.302A>ACC; m.7587T>C"),
      hr(),
      
      downloadButton(ns("download_csv"), "Download Current View", class = "btn-success")
    ),
    
    tags$head(
      tags$style(HTML("
        .dataTables_wrapper { font-size: 13px; }
        table.dataTable { white-space: nowrap; } 
        .card-header { font-weight: 600; font-size: 1.1rem; background-color: #f8f9fa; }
        .dt-buttons { margin-bottom: 10px; }
      "))
    ),
    
    uiOutput(ns("main_content"))
  )
}

allVariantsServer <- function(id, db_path, con) {
  moduleServer(id, function(input, output, session) {
    
    current_view <- reactiveVal("default")
    
    observeEvent(input$view_default, current_view("default"))
    observeEvent(input$view_summary, current_view("summary"))
    observeEvent(input$view_breakup, current_view("breakup"))
    
    observeEvent(db_path(), {
      req(db_path())
      var_path  <- file.path(db_path(), "variants", "*.parquet")
      anno_path <- file.path(db_path(), "annotations", "*.parquet")
      map_path  <- file.path(db_path(), "metadata", "*.parquet")
      
      try({
        dbExecute(con, paste0("CREATE OR REPLACE VIEW variants AS SELECT * FROM read_parquet('", var_path, "')"))
        dbExecute(con, paste0("CREATE OR REPLACE VIEW annotations AS SELECT * FROM read_parquet('", anno_path, "')"))
        dbExecute(con, paste0("CREATE OR REPLACE VIEW sample_map AS SELECT * FROM read_parquet('", map_path, "')"))
      }, silent = TRUE)
    })
    
    base_data <- reactive({
      req(db_path())
      query <- tbl(con, "variants") %>%
        left_join(tbl(con, "annotations"), by = c("VariantID")) %>%
        left_join(tbl(con, "sample_map"), by = c("Sample_ID" = "GRG"))
      
      if (isTruthy(input$search_gene)) {
        search_term <- paste0("%", toupper(input$search_gene), "%")
        query <- query %>% filter(upper(Gene.ensGene) %like% !!search_term)
      }
      
      if (isTruthy(input$search_variant)) {
        search_terms <- str_split(input$search_variant, ";")[[1]] %>% 
          str_trim() %>% 
          toupper() %>% 
          .[. != ""] 
        
        if (length(search_terms) > 0) {
          query <- query %>% filter(upper(VariantID) %in% !!search_terms)
        }
      }
      
      query
    })
    
    data_default <- reactive({
      base_data() %>% 
        mutate(Variant = paste0(CHROM, ".", POS, REF, ">", ALT)) %>%
        dplyr::select(Variant, VariantID, CHROM, POS, REF, ALT, Gene.ensGene, Func.ensGene) %>%
        collect() %>% unique()
    })
    
    data_summary <- reactive({
      raw_data <- base_data() %>% 
        mutate(Variant = paste0(CHROM, ".", POS, REF, ">", ALT)) %>%
        collect()
      
      variant_base <- raw_data %>%
        dplyr::select(VariantID, Variant, CHROM, POS, REF, ALT, Gene.ensGene, Func.ensGene) %>%
        distinct()
      
      # FIX: Added filter(!is.na(GROUP)) before calculating counts
      group_counts <- raw_data %>%
        filter(!is.na(GROUP)) %>% 
        dplyr::select(VariantID, GROUP, NAME) %>%
        distinct() %>% 
        dplyr::group_by(VariantID, GROUP) %>%
        dplyr::summarise(Count = n(), .groups = 'drop') %>%
        tidyr::pivot_wider(names_from = GROUP, values_from = Count, values_fill = 0)
      
      expected_groups <- c("POAG", "Cat", "XFG", "NTG")
      for (col in expected_groups) {
        if (!col %in% names(group_counts)) {
          group_counts[[col]] <- 0
        }
      }
      
      variant_base %>%
        left_join(group_counts, by = "VariantID") %>%
        mutate(across(any_of(expected_groups), ~replace_na(., 0)))
    })
    
    data_breakup <- reactive({
      raw_data <- base_data() %>% 
        mutate(Variant = paste0(CHROM, ".", POS, REF, ">", ALT)) %>%
        collect()
      
      req(nrow(raw_data) > 0)
      raw_data <- raw_data %>% filter(!is.na(NAME))
      
      # FIX: Updated metrics list
      patient_metrics <- c("AD", "AF", "AN", "AC_hom", "AC_het", "AF_hom", "AF_het", "zygosity", "PolyPhen2", "SIFT", "CADD")
      
      for (col in patient_metrics) {
        if (!col %in% names(raw_data)) raw_data[[col]] <- NA
      }
      
      wide_data <- raw_data %>%
        tidyr::pivot_wider(
          id_cols = c("VariantID", "Variant", "Gene.ensGene", "Func.ensGene"),
          names_from = NAME,
          values_from = all_of(patient_metrics),
          names_glue = "{NAME}___{.value}"
        )
      
      patients <- sort(unique(raw_data$NAME))
      col_order <- c("VariantID", "Variant", "Gene.ensGene", "Func.ensGene")
      for (p in patients) {
        col_order <- c(col_order, paste0(p, "___", patient_metrics))
      }
      
      wide_data %>% dplyr::select(all_of(col_order)) %>%
        replace(is.na(.), "Absent")
    })
    
    output$main_content <- renderUI({
      view <- current_view()
      if (view == "default") {
        tagList(h3("All Variants View"), DTOutput(session$ns("table_default")))
        
      } else if (view == "summary") {
        tagList(
          h3("Summary Data"),
          DTOutput(session$ns("table_summary"))
        )
        
      } else if (view == "breakup") {
        tagList(h3("Patient Breakup"), DTOutput(session$ns("table_breakup")))
      }
    })
    
    output$table_default <- renderDT({
      datatable(
        data_default(), 
        class = 'display compact nowrap', 
        width = "100%", 
        options = list(pageLength = 15, scrollX = TRUE, autoWidth = FALSE)
      )
    })
    
    output$table_summary <- renderDT({
      datatable(
        data_summary(), 
        class = 'display compact nowrap', 
        width = "100%", 
        options = list(pageLength = 15, scrollX = TRUE, autoWidth = FALSE)
      )
    })
    
    output$table_breakup <- renderDT({
      wide_data <- data_breakup()
      req(nrow(wide_data) > 0)
      
      # FIX: Updated metrics list in the UI renderer
      patient_metrics <- c("AD", "AF", "AN", "AC_hom", "AC_het", "AF_hom", "AF_het", "zygosity", "PolyPhen2", "SIFT", "CADD")
      all_metric_cols <- names(wide_data)[-(1:4)]
      patients <- unique(sub("___.*", "", all_metric_cols))
      
      th_top <- list(
        tags$th(rowspan = 2, 'VariantID'),
        tags$th(rowspan = 2, 'Variant'),
        tags$th(rowspan = 2, 'Gene.ensGene'),
        tags$th(rowspan = 2, 'Func.ensGene')
      )
      
      for (p in patients) {
        th_top <- c(th_top, list(tags$th(colspan = length(patient_metrics), p, class = 'dt-center')))
      }
      
      th_bottom <- lapply(rep(patient_metrics, length(patients)), tags$th)
      
      sketch <- tags$table(
        class = 'display compact nowrap',
        tags$thead(tags$tr(th_top), tags$tr(th_bottom))
      )
      
      datatable(
        wide_data,
        class = 'display compact nowrap',
        width = "100%",
        container = sketch,
        rownames = FALSE,
        options = list(pageLength = 15, scrollX = TRUE, autoWidth = FALSE)
      )
    })
    
    output$download_csv <- downloadHandler(
      filename = function() { paste0("MitoVar_Export_", current_view(), "_", Sys.Date(), ".csv") },
      content = function(file) {
        view <- current_view()
        export_data <- if (view == "default") {
          data_default()
        } else if (view == "summary") {
          data_summary()
        } else if (view == "breakup") {
          df <- data_breakup()
          names(df) <- gsub("___", "_", names(df))
          df
        }
        write.csv(export_data, file, row.names = FALSE)
      }
    )
    
  })
}