

cl_phv = c("entire_cell_cd19_opal_480_min", "entire_cell_cd19_opal_480_max", 
"entire_cell_cd19_opal_480_std_dev", "entire_cell_cd19_opal_480_total", 
"entire_cell_dapi_dapi_min", "entire_cell_dapi_dapi_max", "entire_cell_dapi_dapi_std_dev", 
"entire_cell_dapi_dapi_total", "entire_cell_autofluorescence_min", 
"entire_cell_autofluorescence_max", "entire_cell_autofluorescence_std_dev", 
"entire_cell_autofluorescence_total")

#' return a vector of colData labels denoting cell phenotypes for EH7312
#' @export
ovVP_phv = function() c("phenotype_cd68", "phenotype_ki67", 
"phenotype_ck", "phenotype_cd19", "phenotype_p_stat3", "phenotype_cd3", 
"phenotype_cd8")

#' return a vector of colData labels denoting cell phenotypes for EH7311
#' @export
luV3_phv = c("phenotype_ck", "phenotype_cd8", "phenotype_cd14", "phenotype_other", 
"phenotype_cd19", "phenotype_cd4")


#' produce a factor with 'positive' cell type values for SpatialExperiment annotated like EH7312
#' @param spe SpatialExperiment instance
#' @param inds numeric() indices of 'phenotype' variables in `colData(spe)`
#' @param opts character() prefixes of values for phenotype values
#' @note This is a speculative implementation.  If multiple factors for a cell
#' are positive, the last one (in the ordering of `opts`) is used.
#' @examples
#' requireNamespace("SummarizedExperiment")
#' data(litov)
#' if ("multicol" %in% names(SummarizedExperiment::colData(litov))) {
#'   wh = which(names(SummarizedExperiment::colData(litov)) == "multicol")
#'   colData(litov) = SummarizedExperiment::colData(litov)[,-wh]
#'   }
#' litov = add_multicol(litov)
#' table(litov$multicol, exclude=NULL)
#' @export
add_multicol = function(spe, inds=188:194, opts = c("CD68", "CK", "CD19", "CD3", "CD8")) {
 mymat = data.matrix(colData(spe)[,inds])
 nopts = paste0(opts, rep(c("+", "-"), each=5))
 nina = function(x) if (length(x)==0) NA else x
 zz = apply(mymat, 1, function(x) last(match(nina(intersect(x, nopts[1:5])), nopts[1:5], nomatch=NA)))
 spe$multicol = factor(opts[zz])
 spe
}

#' use ggplot to display cell phenotype over positions
#' @param spe SpatialExperiment instance
#' @param colvar character(1) phenotype label
#' @param legend_glyph_size numeric(1) passed to `override.aes = list(size=)` in `guides(colour...` defaults to 5
#' @param legend_font_size numeric(1) passed to `theme(legend.text=element_text(size=))`, defaults to 14
#' @param \dots passed to geom_point
#' @examples
#' data(litov)
#' show_core(litov)
#' @export
show_core = function (spe, colvar = "phenotype_cd8", legend_glyph_size = 5, legend_font_size = 14, 
    ...) 
{
    stopifnot(colvar %in% names(colData(spe)))
    xy = spatialCoords(spe)
    datf = data.frame(x = xy[, 1], y = xy[, 2])
    datf$pheno = factor(colData(spe)[, colvar])
    ggplot(datf, aes(x = x, y = y, colour = pheno)) + geom_point(...) + 
        labs(title = spe$sample_id[1]) + guides(colour = guide_legend(override.aes = list(size=legend_glyph_size))) + theme(legend.text = element_text(size = legend_font_size))
}

show_core_old = function(spe, colvar = "phenotype_cd8", legend_font_size=14, ...) {
  stopifnot(colvar %in% names(colData(spe)))
  xy = spatialCoords(spe)
  datf = data.frame(x=xy[,1], y=xy[,2])
  datf$feat = factor(colData(spe)[, colvar])
  ggplot(datf, aes(x=x, y=y, colour=feat)) + geom_point(...) + labs(title=spe$sample_id[1]) +
     theme(legend.text=element_text(size=legend_font_size))
}
 

extract_core = function(id, spe) {
 stopifnot("nid" %in% names(colData(spe)))
 spe[, which(spe$nid==id)]
}

colNplus = function(datf) {
 ans = do.call(cbind, lapply(datf, function(x) length(grep("\\+$", x))))
 names(ans) = names(datf)
 ans
}


get_summaries = function(spe, phv) {
  dat = as.data.frame(colData(spe)[, c("nid", phv) ])
  sdat = split(dat, dat$nid)
  ans = data.frame(do.call(rbind, lapply(sdat, colNplus)))
  ans$nid = names(sdat)
  ans
}

#' Explore ovVP Vectra Polaris data from ExperimentHub EH7312
#' @import SpatialExperiment
#' @import ggplot2
#' @rawNamespace import("shiny", except=c("renderDataTable", "dataTableOutput"))
#' @import DT
#' @import dplyr
#' @param spe the SpatialExperiment from EH7312
#' @param phv character() vector of phenotype labels, from colData
#' @param add_cd4 logical(1) if TRUE, use CD3+ CD8- to enumerate CD4
#' @param add_multicol logical(1) if TRUE, add a multilevel factor for cell type
#' @note  App permits selection of specific tissue slices on tissue microarray, coloring
#' the points at cell coordinates in `spatialCoords`.  The input SpatialExperiment must
#' include a column 'nid' of numerical slice ids in colData.
#' @examples
#' if (interactive()) {
#'   requireNamespace("ExperimentHub")
#'   ovVP = ExperimentHub::ExperimentHub()[["EH7312"]]
#'   ovVP$nid = as.numeric(factor(ovVP$sample_id))
#'   display_by_core(ovVP, ovVP_phv())
#' }
#' @export
display_by_core = function(spe, phv, add_cd4=TRUE, add_multicol=TRUE) {
 stopifnot("nid" %in% names(colData(spe)))
 spe$cd4 = "CD4-"
 if (add_cd4) {
   kp = which(spe$phenotype_cd3=="CD3+")
   spe$cd4[kp] = ifelse(spe$phenotype_cd8[kp]=="CD8-", "CD4+", "CD4-")
   phv = c(phv, "cd4")
   }
 if (add_multicol) phv = c("multicol", phv)
 ui = fluidPage(
  sidebarLayout(
   sidebarPanel(
    selectInput("core", "core", choices=sort(unique(spe$nid))),
    radioButtons("phenotype", "phenotype", choices=phv),
    actionButton("stopBtn", "Stop app"), width=2
    ),
   mainPanel(
    tabsetPanel(
     tabPanel("core", plotlyOutput("coreview", height="550px")),
     tabPanel("phdata", DT::dataTableOutput("datatab")),
     tabPanel("about", helpText("Data described in PMID 34615692.  We use CD3+ and CD8- to label cells CD4+; all others are labeled CD4-."))
     )
    )
   )
  )
 server = function(input, output) {
  output$coreview = renderPlotly({
    newov = extract_core( input$core, spe )
    tmp = newov
    #save(tmp, file="tmp.rda")
    if (add_multicol) {
      newov = add_multicol(newov)
      }
    ggplotly(show_core(newov, input$phenotype))
    })
  output$datatab = DT::renderDataTable({
    if (!exists("ovVP_summaries")) data("ovVP_summaries", package="BiocSpatialDemos")
    ovVP_summaries
    })
  observe({
      if (input$stopBtn > 0) {
        isolate({
          stopApp(returnValue = 0)
        })
      }
    })

 }
 runApp(list(ui=ui, server=server))
}
   
   
    
