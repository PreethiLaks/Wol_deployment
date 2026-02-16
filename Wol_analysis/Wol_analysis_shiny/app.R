#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)

# app.R
# Shiny app: vector suitability heatmaps (terra + sf) with PfPR/ITN filters

library(shiny)
library(terra)
library(sf)

#### ==== 0) EDIT THESE PATHS ==== ####
paths <- list(
  ara       = "C:/Users/LakshminarayananPree/Downloads/2010_Anopheles_arabiensis/2010_Anopheles_arabiensis.tif",
  col       = "C:/Users/LakshminarayananPree/Downloads/2017_Anopheles_coluzzii.Mean_Decompressed (1)/2017_Anopheles_coluzzii.Mean_Decompressed.geotiff",
  mou       = "C:/Users/LakshminarayananPree/Downloads/2010_Anopheles_moucheti/2010_Anopheles_moucheti.tif",
  pfpr_dir  = "C:/Users/LakshminarayananPree/Downloads/2024_GBD2023_Global_PfPR_2000 (2)/",
  itn_dir   = "C:/Users/LakshminarayananPree/Downloads/ITN_2000_use_rate_mean/",
  gadm_gpkg = "C:/Users/LakshminarayananPree/Downloads/gadm_410-levels/gadm_410-levels.gpkg"
)

#### ==== Helpers ==== ####
to01 <- function(r) {  # normalize 0–100 → 0–1 if needed
  mx <- suppressWarnings(as.numeric(global(r, "max", na.rm = TRUE)))
  if (is.finite(mx) && mx > 1) r/100 else r
}
mix_max <- function(a, b) app(c(a, b), function(x) pmax(x[1], x[2], na.rm = TRUE))

build_pfpr_mean <- function(dirpath){
  files <- list.files(dirpath, pattern = "\\.tif(f)?$", full.names = TRUE, ignore.case = TRUE)
  validate(need(length(files)>0, paste0("No PfPR .tif files found in: ", dirpath)))
  mean(rast(files), na.rm = TRUE)
}
build_itn_mean <- function(dirpath){
  files <- list.files(dirpath, pattern = "\\.tif(f)?$", full.names = TRUE, ignore.case = TRUE)
  validate(need(length(files)>0, paste0("No ITN .tif files found in: ", dirpath)))
  mean(rast(files), na.rm = TRUE)
}

plot_with_legend_outside <- function(r, v_adm0, v_adm1, main, pal_name="Viridis"){
  af_ext <- ext(v_adm0)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mar = c(4,4,4,8)) # extra right margin for legend
  plot(r,
       col    = rev(hcl.colors(200, pal_name)),
       zlim   = c(0,1),
       colNA  = "white",
       legend = TRUE,
       axes   = TRUE,
       xlim   = c(xmin(af_ext), xmax(af_ext)),
       ylim   = c(ymin(af_ext), ymax(af_ext)),
       main   = main)
  plot(v_adm1, add=TRUE, lwd=0.2, col=NA, border=adjustcolor("grey85",0.9))
  plot(v_adm0, add=TRUE, lwd=0.8, col=NA, border="black")
}

#### ==== UI ==== ####
ui <- fluidPage(
  titlePanel("Vector suitability heat map (terra + sf)"),
  sidebarLayout(
    sidebarPanel(
      selectInput("species", "Species", c("arabiensis","coluzzii","moucheti"), "arabiensis"),
      selectInput("palette", "Palette", c("Viridis","Inferno","Blues","Greens","Plasma"), "Viridis"),
      numericInput("eps", "Presence epsilon (treat ≤ as absent)", value = 0.001, min = 0, step = 0.001),
      hr(),
      h4("Epi filters"),
      selectInput("pfpr_op", "PfPR operator", c("≥"="ge","<"="lt"), "ge"),
      sliderInput("pfpr_thr", "PfPR threshold (%)", min=0, max=100, value=40, step=5),
      selectInput("itn_op", "ITN operator",  c("≥"="ge","<"="lt"), "ge"),
      sliderInput("itn_thr", "ITN threshold (%)", min=0, max=100, value=80, step=5),
      hr(),
      checkboxInput("show_stats", "Show quick stats", value = TRUE),
      actionButton("run", "Render map")
    ),
    mainPanel(
      plotOutput("map", height = 700),
      verbatimTextOutput("stats")
    )
  )
)

#### ==== Server ==== ####
server <- function(input, output, session){
  
  # Load once
  ara <- reactiveVal(rast(paths$ara))
  colu <- reactiveVal(rast(paths$col))
  mou  <- reactiveVal(rast(paths$mou))
  pfpr_mean <- reactiveVal(build_pfpr_mean(paths$pfpr_dir))
  itn_mean  <- reactiveVal(build_itn_mean(paths$itn_dir))
  
  # Borders (GADM v4.1) as sf → spatVector
  admin0 <- reactiveVal( vect(st_read(paths$gadm_gpkg, layer="ADM_0", quiet = TRUE)) )
  admin1 <- reactiveVal( vect(st_read(paths$gadm_gpkg, layer="ADM_1", quiet = TRUE)) )
  
  # Compute suitability for selected species after button press
  suit_r <- eventReactive(input$run, {
    sp <- input$species
    eps <- input$eps
    
    # Select base species raster + align "others" to it
    if (sp == "arabiensis") {
      base <- ara()
      other1 <- project(colu(), base, method="bilinear")
      other2 <- project(mou(),  base, method="bilinear")
    } else if (sp == "coluzzii") {
      base <- colu()
      other1 <- project(ara(), base, method="bilinear")
      other2 <- project(mou(), base, method="bilinear")
    } else {
      base <- mou()
      other1 <- project(ara(), base, method="bilinear")
      other2 <- project(colu(), base, method="bilinear")
    }
    
    # Normalize to 0–1
    base01   <- to01(base)
    other101 <- to01(other1)
    other201 <- to01(other2)
    
    # Mixture = max of other vectors
    other101[is.na(other101)] <- 0
    other201[is.na(other201)] <- 0
    mix  <- mix_max(other101, other201)              # [0,1]
    suit <- app(mix, function(x) 1 - pmin(x, 0.99))  # suitability [0,1], exclusive darkest
    
    # Strict presence mask for base species: base > eps => keep; else NA
    pres_mask <- ifel(base01 > eps, 1, NA)
    suit <- suit * pres_mask  # absence/no-data -> NA
    
    # Epi filters (PfPR & ITN) on the same grid
    pfpr_b <- project(pfpr_mean(), base, method="bilinear")
    itn_b  <- project(itn_mean(),  base, method="bilinear")
    pfpr01 <- to01(pfpr_b)
    itn01  <- to01(itn_b)
    
    pfpr_thr <- input$pfpr_thr / 100
    itn_thr  <- input$itn_thr  / 100
    pf_ok <- if (input$pfpr_op=="ge") pfpr01 >= pfpr_thr else pfpr01 < pfpr_thr
    it_ok <- if (input$itn_op =="ge") itn01  >= itn_thr  else itn01  < itn_thr
    epi_mask <- ifel(pf_ok & it_ok, 1, 0)
    
    suit_epi <- mask(suit, epi_mask, maskvalues = 0)
    
    list(suit = suit_epi, admin0 = admin0(), admin1 = admin1())
  }, ignoreInit = TRUE)
  
  # Map
  output$map <- renderPlot({
    req(suit_r())
    x <- suit_r()
    # Reproject borders to match raster CRS (safest)
    v0 <- if (!identical(crs(x$admin0), crs(x$suit))) project(x$admin0, crs(x$suit)) else x$admin0
    v1 <- if (!identical(crs(x$admin1), crs(x$suit))) project(x$admin1, crs(x$suit)) else x$admin1
    
    pf_sym <- if (input$pfpr_op == "ge") "≥" else "<"
    it_sym <- if (input$itn_op  == "ge") "≥" else "<"
    
    ttl <- sprintf(
      "%s suitability (PfPR %s %d%%, ITN %s %d%%)",
      tools::toTitleCase(input$species), pf_sym, input$pfpr_thr, it_sym, input$itn_thr
    )
    
    plot_with_legend_outside(x$suit, v0, v1, main = ttl, pal_name = input$palette)
  })
  
  # Stats (fixed: extract numerics from terra::global)
  output$stats <- renderPrint({
    req(input$run, suit_r(), input$show_stats)
    s <- suit_r()$suit
    n_cells <- as.numeric(terra::global(!is.na(s), "sum",    na.rm = TRUE))
    mean_v  <- as.numeric(terra::global(s,          "mean",  na.rm = TRUE))
    med_v   <- as.numeric(terra::global(s,          "median",na.rm = TRUE))
    
    cat("Cells (non-NA):", n_cells, "\n")
    cat("Mean suitability:", round(mean_v, 3), "\n")
    cat("Median suitability:", round(med_v, 3), "\n")
  })
}

shinyApp(ui, server)

