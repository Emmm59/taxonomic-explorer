# Interactive Hierarchical Taxonomic Explorer
#===============================================================================

# Load required packages
suppressPackageStartupMessages({
  library(shiny)
  library(tidyverse)
  library(plotly)
  library(glue)
  library(DT)
  library(viridis)
  library(htmlwidgets)
  library(webshot)
  library(readxl)
})

# Set global parameters
set.seed(42)
font_family <- "Helvetica Neue"
plot_width <- 1200
plot_height <- 800

distribution_colors <- c(
  "Nativa" = "#ffb84d",
  "Introducida" = "#E54F6D",
  "Endémica" = "#85C1A5"
)

allelopathy_colors <- c(
  "Sí" = "#EDF8FB",
  "No" = "#EDF8FB",
  "No Reportado" = "#7f7f7f"
)

#-------------------------------------------------------------------------------
# Data Processing Functions
#-------------------------------------------------------------------------------

validate_data_structure <- function(data) {
  required_cols <- c("Family", "Genus", "Species", "mxc_status", "Alelopatía", "lifeform_description", "climate_description")
  missing_cols <- setdiff(required_cols, names(data))
  
  if (length(missing_cols) > 0) {
    stop("Input data missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Validate critical fields
  if (any(is.na(data$Family) | data$Family == "")) {
    warning("Removing records with missing family names")
    data <- data %>% filter(!is.na(Family) & Family != "")
  }
  
  if (any(is.na(data$Species) | data$Species == "")) {
    stop("Species names cannot be missing - please fix data source")
  }
  
  return(data)
}

clean_distribution_data <- function(data) {
  data %>%
    mutate(
      Genus = str_trim(Genus),
      Species = str_trim(Species),
      scientific_name = paste(Genus, Species),
      lifeform_description = str_to_sentence(lifeform_description),
      mxc_status = case_when(
        str_to_lower(mxc_status) %in% c("nativa", "native") ~ "Nativa",
        str_to_lower(mxc_status) %in% c("introducida", "introduced") ~ "Introducida",
        str_to_lower(mxc_status) %in% c("endémica", "endemic", "endemica") ~ "Endémica"
      ),
      Alelopatía = case_when(
        str_to_lower(Alelopatía) %in% c("si", "sí", "yes", "y", "1") ~ "Sí",
        str_to_lower(Alelopatía) %in% c("no", "n", "0") ~ "No",
        TRUE ~ "No Reportado"
      )
    ) %>%
    distinct(
      Family, Genus, scientific_name, mxc_status, Alelopatía,
      lifeform_description, climate_description
    )
}

calculate_totals <- function(data) {
  data %>%
    summarize(
      total_species = n_distinct(scientific_name),
      total_native = sum(mxc_status == "Nativa", na.rm = TRUE),
      total_introduced = sum(mxc_status == "Introducida", na.rm = TRUE),
      total_endemic = sum(mxc_status == "Endémica", na.rm = TRUE),
      total_allelopathy_yes = sum(Alelopatía == "Sí", na.rm = TRUE),
      total_allelopathy_no = sum(Alelopatía == "No", na.rm = TRUE)
    ) %>%
    as.list()
}

prepare_main_hierarchy <- function(data, hierarchy_type) {
  # Return empty dataframe if no data
  if (nrow(data) == 0) {
    return(data.frame())
  }
  
  if (hierarchy_type == "Taxonómica (Familia > Género)") {
    # Family > Genus > Species hierarchy
    n_families <- n_distinct(data$Family)
    family_colors <- viridis(n_families, option = "D", direction = 1)
    
    family_data <- data %>%
      group_by(Family) %>%
      summarize(
        total_species = n_distinct(scientific_name),
        native = sum(mxc_status == "Nativa"),
        introduced = sum(mxc_status == "Introducida"),
        endemic = sum(mxc_status == "Endémica"),
        .groups = "drop"
      ) %>%
      mutate(
        id = paste0("family_", Family),
        parent = "",
        label = glue("{Family} ({native}|{introduced}|{endemic})"),
        value = total_species,
        color = family_colors[rank(-total_species, ties.method = "first")],
        level = "family"
      )
    
    genus_data <- data %>%
      group_by(Family, Genus) %>%
      summarize(
        total_species = n_distinct(scientific_name),
        .groups = "drop"
      ) %>%
      mutate(
        id = paste0("genus_", Family, "_", Genus),
        parent = paste0("family_", Family),
        label = Genus,
        value = total_species,
        level = "genus"
      ) %>%
      group_by(Family) %>%
      mutate(
        color = colorRampPalette(c("#5E4FA2", "#3288BD"))(n())[rank(-total_species)]
      ) %>%
      ungroup()
    
    species_data <- data %>%
      mutate(
        id = paste0("species_", str_replace_all(scientific_name, " ", "_")),
        parent = paste0("genus_", Family, "_", Genus),
        label = scientific_name,
        value = 1,
        color = distribution_colors[mxc_status],
        level = "species"
      )
    
    hierarchical_data <- bind_rows(
      family_data,
      genus_data,
      species_data
    )
  } 
  else if (hierarchy_type == "Formas de Vida") {
    # Lifeform-based hierarchy: Lifeform → Family → Genus → Species
    n_lifeforms <- n_distinct(data$lifeform_description)
    lifeform_colors <- viridis(n_lifeforms, option = "D", direction = 1)
    
    lifeform_summary <- data %>%
      group_by(lifeform_description) %>%
      summarise(
        total_species = n_distinct(scientific_name),
        native = sum(mxc_status == "Nativa"),
        introduced = sum(mxc_status == "Introducida"),
        endemic = sum(mxc_status == "Endémica"),
        .groups = "drop"
      ) %>%
      mutate(
        id = paste0("lifeform_", lifeform_description),
        parent = "",
        label = glue("{lifeform_description} ({native}|{introduced}|{endemic})"),
        value = total_species,
        color = lifeform_colors[rank(-total_species, ties.method = "first")],
        level = "lifeform"
      )
    
    family_summary <- data %>%
      group_by(lifeform_description, Family) %>%
      summarise(
        total_species = n_distinct(scientific_name),
        native = sum(mxc_status == "Nativa"),
        introduced = sum(mxc_status == "Introducida"),
        endemic = sum(mxc_status == "Endémica"),
        .groups = "drop"
      ) %>%
      mutate(
        id = paste0("family_", Family, "_in_", lifeform_description),
        parent = paste0("lifeform_", lifeform_description),
        label = glue("{Family} ({native}|{introduced}|{endemic})"),
        value = total_species,
        level = "family"
      ) %>%
      group_by(lifeform_description) %>%
      mutate(
        color = colorRampPalette(c("#5E4FA2", "#3288BD"))(n())[rank(-total_species)]
      ) %>%
      ungroup()
    
    # Genus level
    genus_summary <- data %>%
      group_by(lifeform_description, Family, Genus) %>%
      summarise(
        total_species = n_distinct(scientific_name),
        .groups = "drop"
      ) %>%
      mutate(
        id = paste0("genus_", Genus, "_in_", Family, "_", lifeform_description),
        parent = paste0("family_", Family, "_in_", lifeform_description),
        label = Genus,
        value = total_species,
        level = "genus"
      ) %>%
      group_by(lifeform_description, Family) %>%
      mutate(
        color = colorRampPalette(c("#66C2A5", "#3288BD"))(n())[rank(-total_species)]
      ) %>%
      ungroup()
    
    species_data <- data %>%
      mutate(
        id = paste0("species_", str_replace_all(scientific_name, " ", "_")),
        parent = paste0("genus_", Genus, "_in_", Family, "_", lifeform_description),
        label = scientific_name,
        value = 1,
        color = distribution_colors[mxc_status],
        level = "species"
      )
    
    hierarchical_data <- bind_rows(
      lifeform_summary,
      family_summary,
      genus_summary,
      species_data
    )
  }
  else { # Climate-based hierarchy: Climate → Family → Genus → Species
    n_climates <- n_distinct(data$climate_description)
    climate_colors <- viridis(n_climates, option = "D", direction = 1)
    
    climate_summary <- data %>%
      group_by(climate_description) %>%
      summarise(
        total_species = n_distinct(scientific_name),
        native = sum(mxc_status == "Nativa"),
        introduced = sum(mxc_status == "Introducida"),
        endemic = sum(mxc_status == "Endémica"),
        .groups = "drop"
      ) %>%
      mutate(
        id = paste0("climate_", climate_description),
        parent = "",
        label = glue("{climate_description} ({native}|{introduced}|{endemic})"),
        value = total_species,
        color = climate_colors[rank(-total_species, ties.method = "first")],
        level = "climate"
      )
    
    family_summary <- data %>%
      group_by(climate_description, Family) %>%
      summarise(
        total_species = n_distinct(scientific_name),
        native = sum(mxc_status == "Nativa"),
        introduced = sum(mxc_status == "Introducida"),
        endemic = sum(mxc_status == "Endémica"),
        .groups = "drop"
      ) %>%
      mutate(
        id = paste0("family_", Family, "_in_", climate_description),
        parent = paste0("climate_", climate_description),
        label = glue("{Family} ({native}|{introduced}|{endemic})"),
        value = total_species,
        level = "family"
      ) %>%
      group_by(climate_description) %>%
      mutate(
        color = colorRampPalette(c("#5E4FA2", "#3288BD"))(n())[rank(-total_species)]
      ) %>%
      ungroup()
    
    # Genus level
    genus_summary <- data %>%
      group_by(climate_description, Family, Genus) %>%
      summarise(
        total_species = n_distinct(scientific_name),
        .groups = "drop"
      ) %>%
      mutate(
        id = paste0("genus_", Genus, "_in_", Family, "_", climate_description),
        parent = paste0("family_", Family, "_in_", climate_description),
        label = Genus,
        value = total_species,
        level = "genus"
      ) %>%
      group_by(climate_description, Family) %>%
      mutate(
        color = colorRampPalette(c("#66C2A5", "#3288BD"))(n())[rank(-total_species)]
      ) %>%
      ungroup()
    
    species_data <- data %>%
      mutate(
        id = paste0("species_", str_replace_all(scientific_name, " ", "_")),
        parent = paste0("genus_", Genus, "_in_", Family, "_", climate_description),
        label = scientific_name,
        value = 1,
        color = distribution_colors[mxc_status],
        level = "species"
      )
    
    hierarchical_data <- bind_rows(
      climate_summary,
      family_summary,
      genus_summary,
      species_data
    )
  }
  
  # Add tooltips
  hierarchical_data %>%
    mutate(
      tooltip = case_when(
        level == "family" ~ glue(
          "<b>{Family}</b><br>",
          "Total especies: {value}<br>",
          "<span style='color:{distribution_colors['Nativa']}'>Nativas: {native}</span><br>",
          "<span style='color:{distribution_colors['Introducida']}'>Introducidas: {introduced}</span><br>",
          "<span style='color:{distribution_colors['Endémica']}'>Endémicas: {endemic}</span>"
        ),
        level == "genus" ~ glue(
          "<b>{label}</b><br>",
          "<span style='color:#888888'>Total especies: {value}</span>"
        ),
        level == "lifeform" ~ glue(
          "<b>{lifeform_description}</b><br>",
          "Total especies: {value}<br>",
          "<span style='color:{distribution_colors['Nativa']}'>Nativas: {native}</span><br>",
          "<span style='color:{distribution_colors['Introducida']}'>Introducidas: {introduced}</span><br>",
          "<span style='color:{distribution_colors['Endémica']}'>Endémicas: {endemic}</span>"
        ),
        level == "climate" ~ glue(
          "<b>{climate_description}</b><br>",
          "Total especies: {value}<br>",
          "<span style='color:{distribution_colors['Nativa']}'>Nativas: {native}</span><br>",
          "<span style='color:{distribution_colors['Introducida']}'>Introducidas: {introduced}</span><br>",
          "<span style='color:{distribution_colors['Endémica']}'>Endémicas: {endemic}</span>"
        ),
        level == "species" ~ glue(
          "<span style='color:{distribution_colors[mxc_status]}'><b>{label}</b></span><br>",
          "Distribución: {mxc_status}<br>",
          "Alelopatía: {Alelopatía}<br>",
          "Forma de vida: {lifeform_description}<br>",
          "Clima: {climate_description}"
        )
      )
    )
}

create_treemap <- function(hierarchical_data, totals, title_suffix = "") {
  # Calculate allelopathy stats for title
  allelopathy_stats <- ""
  if ("Alelopatía" %in% names(hierarchical_data)) {
    allelopathy_yes <- sum(hierarchical_data$Alelopatía == "Sí", na.rm = TRUE)
    allelopathy_no <- sum(hierarchical_data$Alelopatía == "No", na.rm = TRUE)
    allelopathy_stats <- glue(
      " | Con alelopatía: {allelopathy_yes}",
      " | Sin alelopatía: {allelopathy_no}"
    )
  }
  
  title_text <- glue(
    "<b>DISTRIBUCIÓN TAXONÓMICA DE ESPECIES VEGETALES</b>  | ",
    "<span style='font-size:14px; color:#555555; font-weight:normal; margin-left:20px'>",
    "Total especies: {totals$total_species} | ",
    "Nativas: {totals$total_native} | ",
    "Introducidas: {totals$total_introduced} | ",
    "Endémicas: {totals$total_endemic}",
    "{allelopathy_stats}</span>"
  )
  
  plot_ly(
    data = hierarchical_data,
    type = "treemap",
    ids = ~id,
    labels = ~label,
    parents = ~parent,
    values = ~value,
    text = ~tooltip,
    hoverinfo = "text",
    marker = list(colors = ~color),
    branchvalues = "total",
    pathbar = list(visible = TRUE, thickness = 25),
    textfont = list(family = font_family, size = 14),
    insidetextfont = list(family = font_family, size = 40, color = "white"),
    outsidetextfont = list(family = font_family, size = 12, color = "#2c3e50")
  ) %>%
    layout(
      title = list(
        text = title_text,
        x = 0.02,
        y = 0.98,
        font = list(family = font_family, size = 22, color = "#2c3e50"),
        xanchor = "left"
      ),
      margin = list(t = 73, l = 10, r = 10, b = 10),
      paper_bgcolor = "#f8f9fa",
      hoverlabel = list(
        bgcolor = "rgba(255,255,255,0.9)",
        bordercolor = "#cccccc",
        font = list(family = font_family, size = 12, color = "#333333")
      ),
      annotations = list(
        list(
          x = 0, y = -0.1,
          text = "Haga clic en una categoría para ver niveles inferiores",
          showarrow = FALSE,
          xref = "paper",
          yref = "paper",
          font = list(size = 11, color = "#6c757d")
        )
      )
    ) %>%
    config(
      displayModeBar = TRUE,
      modeBarButtonsToRemove = c("zoom2d", "select2d", "lasso2d", "autoScale2d"),
      toImageButtonOptions = list(
        format = "png",
        filename = "taxonomic_distribution",
        width = plot_width,
        height = plot_height
      )
    )
}

#-------------------------------------------------------------------------------
# Shiny Application
#-------------------------------------------------------------------------------

# Load and preprocess data
raw_data <- readxl::read_excel("BaseTesis.xlsx")
validated_data <- validate_data_structure(raw_data)
clean_data <- clean_distribution_data(validated_data)

# Define UI
ui <- fluidPage(
  titlePanel("Explorador Taxonómico de Especies Vegetales"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Filtros de Datos"),
      checkboxGroupInput(
        "status_filter", "Distribución:",
        choices = c("Nativa", "Introducida", "Endémica"),
        selected = c("Nativa", "Introducida", "Endémica")
      ),
      selectizeInput(
        "family_filter", "Familias:",
        choices = sort(unique(clean_data$Family)),
        multiple = TRUE,
        options = list(placeholder = 'Todas las familias')
      ),
      selectizeInput(
        "lifeform_filter", "Formas de Vida:",
        choices = sort(unique(clean_data$lifeform_description)),
        multiple = TRUE,
        options = list(placeholder = 'Todas las formas de vida')
      ),
      selectizeInput(
        "climate_filter", "Tipos de Clima:",
        choices = sort(unique(clean_data$climate_description)),
        multiple = TRUE,
        options = list(placeholder = 'Todos los climas')
      ),
      radioButtons(
        "allelopathy_filter", "Allelopathy:",
        choices = c("Todas", "Sí", "No", "No Reportado"),
        selected = "Todas"
      ),
      radioButtons(
        "hierarchy_type", "Tipo de Jerarquía:",
        choices = c("Taxonómica (Familia > Género)", 
                    "Formas de Vida", 
                    "Climas"),
        selected = "Taxonómica (Familia > Género)"
      ),
      actionButton("apply_filters", "Aplicar Filtros", class = "btn-primary"),
      hr(),
      downloadButton("download_plot", "Descargar Visualización")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel(
          "Visualización Taxonómica",
          plotlyOutput("main_treemap", height = "800px")
        ),
        tabPanel(
          "Datos",
          DTOutput("data_table")
        )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  # Reactive data processing
  filtered_data <- eventReactive(input$apply_filters, {
    data <- clean_data
    
    # Apply status filters
    if (!is.null(input$status_filter)) {
      data <- data %>% filter(mxc_status %in% input$status_filter)
    }
    
    # Apply family filters
    if (!is.null(input$family_filter)) {
      data <- data %>% filter(Family %in% input$family_filter)
    }
    
    # Apply lifeform filters
    if (!is.null(input$lifeform_filter)) {
      data <- data %>% filter(lifeform_description %in% input$lifeform_filter)
    }
    
    # Apply climate filters
    if (!is.null(input$climate_filter)) {
      data <- data %>% filter(climate_description %in% input$climate_filter)
    }
    
    # Apply allelopathy filter
    if (input$allelopathy_filter != "Todas") {
      data <- data %>% filter(Alelopatía == input$allelopathy_filter)
    }
    
    data
  })
  
  # Update family choices based on other filters
  observe({
    # Get current filters
    status <- if (!is.null(input$status_filter)) input$status_filter else unique(clean_data$mxc_status)
    lifeforms <- if (!is.null(input$lifeform_filter)) input$lifeform_filter else unique(clean_data$lifeform_description)
    climates <- if (!is.null(input$climate_filter)) input$climate_filter else unique(clean_data$climate_description)
    
    # Filter data based on current selections
    filtered_families <- clean_data %>%
      filter(mxc_status %in% status,
             lifeform_description %in% lifeforms,
             climate_description %in% climates) %>%
      pull(Family) %>%
      unique() %>%
      sort()
    
    # Update family choices
    updateSelectizeInput(
      session, 
      "family_filter",
      choices = filtered_families,
      selected = input$family_filter
    )
  })
  
  # Reactive plot
  main_plot <- reactive({
    data <- filtered_data()
    totals <- calculate_totals(data)
    
    hierarchical_data <- prepare_main_hierarchy(data, input$hierarchy_type)
    create_treemap(hierarchical_data, totals)
  })
  
  # UI outputs
  output$main_treemap <- renderPlotly({
    plot <- main_plot()
    
    # Handle empty dataset
    if (nrow(filtered_data()) == 0) {
      return(plotly_empty(type = "scatter") %>%
               layout(
                 title = "No hay datos disponibles con los filtros actuales",
                 plot_bgcolor = "#f8f9fa"
               ))
    }
    plot
  })
  
  output$data_table <- renderDT({
    data <- filtered_data()
    
    # Handle empty dataset
    if (nrow(data) == 0) {
      return(data.frame(Mensaje = "No hay datos disponibles con los filtros actuales"))
    }
    
    data %>%
      select(
        Familia = Family,
        Género = Genus,
        Especie = scientific_name,
        Distribución = mxc_status,
        Alelopatía,
        `Forma de vida` = lifeform_description,
        Clima = climate_description
      )
  }, options = list(pageLength = 10, scrollX = TRUE))
  
  # Download handler
  output$download_plot <- downloadHandler(
    filename = function() {
      "taxonomic_distribution.png"
    },
    content = function(file) {
      showNotification("Generando imagen...", duration = NULL, id = "download_msg")
      on.exit(removeNotification("download_msg"))
      
      plot <- main_plot()
      
      # Handle empty dataset
      if (nrow(filtered_data()) == 0) {
        stop("No hay datos disponibles para descargar")
      }
      
      temp_html <- tempfile(fileext = ".html")
      temp_png <- tempfile(fileext = ".png")
      
      htmlwidgets::saveWidget(
        widget = as_widget(plot), 
        file = temp_html, 
        selfcontained = TRUE
      )
      webshot::webshot(
        url = temp_html, 
        file = temp_png,
        vwidth = plot_width, 
        vheight = plot_height, 
        delay = 3
      )
      file.copy(temp_png, file, overwrite = TRUE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)