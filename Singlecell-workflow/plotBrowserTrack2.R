plotBrowserTrack2 <- function (ArchRProj = NULL, region = NULL, groupBy = "Clusters", 
          useGroups = NULL, plotSummary = c("bulkTrack", "featureTrack", 
                                            "loopTrack", "geneTrack"), sizes = c(10, 1.5, 3, 4), 
          features = getPeakSet(ArchRProj), loops = getCoAccessibility(ArchRProj), 
          geneSymbol = NULL, useMatrix = NULL, log2Norm = TRUE, upstream = 50000, 
          downstream = 50000, tileSize = 250, minCells = 25, normMethod = "ReadsInTSS", 
          threads = getArchRThreads(), ylim = NULL, pal = NULL, baseSize = 7, 
          scTileSize = 0.5, scCellsMax = 100, borderWidth = 0.4, tickWidth = 0.4, 
          facetbaseSize = 7, geneAnnotation = getGeneAnnotation(ArchRProj), 
          title = "", verbose = TRUE, logFile = createLogFile("plotBrowserTrack")) 
{
    ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProj")
    ArchR:::.validInput(input = region, name = "region", valid = c("granges", 
                                                           "null"))
    ArchR:::.validInput(input = groupBy, name = "groupBy", valid = "character")
    ArchR:::.validInput(input = useGroups, name = "useGroups", valid = c("character", 
                                                                 "null"))
    ArchR:::.validInput(input = plotSummary, name = "plotSummary", valid = "character")
    ArchR:::.validInput(input = sizes, name = "sizes", valid = "numeric")
    ArchR:::.validInput(input = features, name = "features", valid = c("granges", 
                                                               "grangeslist", "null"))
    ArchR:::.validInput(input = loops, name = "loops", valid = c("granges", 
                                                         "grangeslist", "null"))
    ArchR:::.validInput(input = geneSymbol, name = "geneSymbol", valid = c("character", 
                                                                   "null"))
    ArchR:::.validInput(input = useMatrix, name = "useMatrix", valid = c("character", 
                                                                 "null"))
    ArchR:::.validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
    ArchR:::.validInput(input = upstream, name = "upstream", valid = c("integer"))
    ArchR:::.validInput(input = downstream, name = "downstream", valid = c("integer"))
    ArchR:::.validInput(input = tileSize, name = "tileSize", valid = c("integer"))
    ArchR:::.validInput(input = minCells, name = "minCells", valid = c("integer"))
    ArchR:::.validInput(input = normMethod, name = "normMethod", valid = c("character"))
    ArchR:::.validInput(input = threads, name = "threads", valid = c("integer"))
    ArchR:::.validInput(input = ylim, name = "ylim", valid = c("numeric", 
                                                       "null"))
    ArchR:::.validInput(input = pal, name = "pal", valid = c("palette", 
                                                     "null"))
    ArchR:::.validInput(input = baseSize, name = "baseSize", valid = "numeric")
    ArchR:::.validInput(input = scTileSize, name = "scTileSize", valid = "numeric")
    ArchR:::.validInput(input = scCellsMax, name = "scCellsMax", valid = "integer")
    ArchR:::.validInput(input = borderWidth, name = "borderWidth", valid = "numeric")
    ArchR:::.validInput(input = tickWidth, name = "tickWidth", valid = "numeric")
    ArchR:::.validInput(input = facetbaseSize, name = "facetbaseSize", 
                valid = "numeric")
    geneAnnotation <- ArchR:::.validGeneAnnotation(geneAnnotation)
    ArchR:::.validInput(input = title, name = "title", valid = "character")
    tstart <- Sys.time()
    ArchR:::.startLogging(logFile = logFile)
    ArchR:::.logThis(mget(names(formals()), sys.frame(sys.nframe())), 
             "plotBrowserTrack Input-Parameters", logFile = logFile)
    ArchR:::.logDiffTime("Validating Region", t1 = tstart, verbose = verbose, 
                 logFile = logFile)
    if (is.null(region)) {
        if (!is.null(geneSymbol)) {
            region <- geneAnnotation$genes
            region <- region[which(tolower(mcols(region)$symbol) %in% 
                                       tolower(geneSymbol))]
            region <- region[order(match(tolower(mcols(region)$symbol), 
                                         tolower(geneSymbol)))]
            print(region)
            region <- GenomicRanges::resize(region, 1, "start")
            strand(region) <- "*"
            region <- extendGR(region, upstream = upstream, 
                               downstream = downstream)
        }
    }
    region <- ArchR:::.validGRanges(region)
    ArchR:::.logThis(region, "region", logFile = logFile)
    if (is.null(geneSymbol)) {
        useMatrix <- NULL
    }
    if (!is.null(useMatrix)) {
        featureMat <- ArchR:::.getMatrixValues(ArchRProj = ArchRProj, 
                                       matrixName = useMatrix, name = mcols(region)$symbol)
        if (log2Norm) {
            featureMat <- log2(featureMat + 1)
        }
        featureMat <- data.frame(t(featureMat))
        featureMat$Group <- getCellColData(ArchRProj, groupBy, 
                                           drop = FALSE)[rownames(featureMat), 1]
    }
    ggList <- lapply(seq_along(region), function(x) {
        plotList <- list()
        if ("bulktrack" %in% tolower(plotSummary)) {
            ArchR:::.logDiffTime(sprintf("Adding Bulk Tracks (%s of %s)", 
                                 x, length(region)), t1 = tstart, verbose = verbose, 
                         logFile = logFile)
            plotList$bulktrack <- .bulkTracks(ArchRProj = ArchRProj, 
                                              region = region[x], tileSize = tileSize, groupBy = groupBy, 
                                              threads = threads, minCells = minCells, pal = pal, 
                                              ylim = ylim, baseSize = baseSize, borderWidth = borderWidth, 
                                              tickWidth = tickWidth, facetbaseSize = facetbaseSize, 
                                              normMethod = normMethod, geneAnnotation = geneAnnotation, 
                                              title = title, useGroups = useGroups, tstart = tstart, 
                                              logFile = logFile) + theme(plot.margin = unit(c(0.35, 
                                                                                              0.75, 0.35, 0.75), "cm"))
        }
        if ("sctrack" %in% tolower(plotSummary)) {
            ArchR:::.logDiffTime(sprintf("Adding SC Tracks (%s of %s)", 
                                 x, length(region)), t1 = tstart, verbose = verbose, 
                         logFile = logFile)
            plotList$sctrack <- ArchR:::.scTracks(ArchRProj = ArchRProj, 
                                          region = region[x], tileSize = tileSize, groupBy = groupBy, 
                                          threads = threads, minCells = 5, maxCells = scCellsMax, 
                                          pal = pal, baseSize = baseSize, borderWidth = borderWidth, 
                                          tickWidth = tickWidth, scTileSize = scTileSize, 
                                          facetbaseSize = facetbaseSize, geneAnnotation = geneAnnotation, 
                                          title = title, useGroups = useGroups, tstart = tstart, 
                                          logFile = logFile) + theme(plot.margin = unit(c(0.35, 
                                                                                          0.75, 0.35, 0.75), "cm"))
        }
        if ("featuretrack" %in% tolower(plotSummary)) {
            if (!is.null(features)) {
                ArchR:::.logDiffTime(sprintf("Adding Feature Tracks (%s of %s)", 
                                     x, length(region)), t1 = tstart, verbose = verbose, 
                             logFile = logFile)
                plotList$featuretrack <- ArchR:::.featureTracks(features = features, 
                                                        region = region[x], facetbaseSize = facetbaseSize, 
                                                        hideX = TRUE, title = "Peaks", logFile = logFile) + 
                    theme(plot.margin = unit(c(0.1, 0.75, 0.1, 
                                               0.75), "cm"))
            }
        }
        if ("looptrack" %in% tolower(plotSummary)) {
            if (!is.null(loops)) {
                ArchR:::.logDiffTime(sprintf("Adding Loop Tracks (%s of %s)", 
                                     x, length(region)), t1 = tstart, verbose = verbose, 
                             logFile = logFile)
                plotList$looptrack <- ArchR:::.loopTracks(loops = loops, 
                                                  region = region[x], facetbaseSize = facetbaseSize, 
                                                  hideX = TRUE, hideY = TRUE, title = "Loops", 
                                                  logFile = logFile) + theme(plot.margin = unit(c(0.1, 
                                                                                                  0.75, 0.1, 0.75), "cm"))
            }
        }
        if ("genetrack" %in% tolower(plotSummary)) {
            ArchR:::.logDiffTime(sprintf("Adding Gene Tracks (%s of %s)", 
                                 x, length(region)), t1 = tstart, verbose = verbose, 
                         logFile = logFile)
            plotList$genetrack <- ArchR:::.geneTracks(geneAnnotation = geneAnnotation, 
                                              region = region[x], facetbaseSize = facetbaseSize, 
                                              title = "Genes", logFile = logFile) + theme(plot.margin = unit(c(0.1, 
                                                                                                               0.75, 0.1, 0.75), "cm"))
        }
        plotSummary <- tolower(plotSummary)
        names(sizes) <- plotSummary
        sizes <- sizes[order(plotSummary)]
        plotSummary <- plotSummary[order(plotSummary)]
        sizes <- sizes[tolower(names(plotList))]
        if (!is.null(useMatrix)) {
            suppressWarnings(ArchR:::.combinedFeaturePlot(plotList = plotList, 
                                                  log2Norm = log2Norm, featureMat = featureMat, 
                                                  feature = region[x]$symbol[[1]], useMatrix = useMatrix, 
                                                  pal = pal, sizes = sizes, baseSize = baseSize, 
                                                  facetbaseSize = facetbaseSize, borderWidth = borderWidth, 
                                                  tickWidth = tickWidth))
        }
        else {
            ArchR:::.logThis(names(plotList), sprintf("(%s of %s) names(plotList)", 
                                              x, length(region)), logFile = logFile)
            ArchR:::.logThis(sizes, sprintf("(%s of %s) sizes", x, length(region)), 
                     logFile = logFile)
            ArchR:::.logDiffTime("Plotting", t1 = tstart, verbose = verbose, 
                         logFile = logFile)
            tryCatch({
                suppressWarnings(ggAlignPlots(plotList = plotList, 
                                              sizes = sizes, draw = FALSE))
            }, error = function(e) {
                ArchR:::.logMessage("Error with plotting, diagnosing each element", 
                            verbose = TRUE, logFile = logFile)
                for (i in seq_along(plotList)) {
                    tryCatch({
                        print(plotList[[i]])
                    }, error = function(f) {
                        .logError(f, fn = names(plotList)[i], info = "", 
                                  errorList = NULL, logFile = logFile)
                    })
                }
                ArchR:::.logError(e, fn = "ggAlignPlots", info = "", 
                          errorList = NULL, logFile = logFile)
            })
        }
    })
    if (!is.null(mcols(region)$symbol)) {
        names(ggList) <- mcols(region)$symbol
    }
    else {
        if (length(ggList) == 1) {
            ggList <- ggList[[1]]
        }
    }
    ArchR:::.endLogging(logFile = logFile)
    ggList
}



.bulkTracks <- function (ArchRProj = NULL, region = NULL, tileSize = 100, minCells = 25, 
          groupBy = "Clusters", useGroups = NULL, normMethod = "ReadsInTSS", 
          threads = 1, ylim = NULL, baseSize = 7, borderWidth = 0.4, 
          tickWidth = 0.4, facetbaseSize = 7, geneAnnotation = getGeneAnnotation(ArchRProj), 
          title = "", pal = NULL, tstart = NULL, verbose = FALSE, 
          logFile = NULL) 
{
    if (is.null(tstart)) {
        tstart <- Sys.time()
    }
    df <-
        ArchR:::.groupRegionSumArrows(
            ArchRProj = ArchRProj,
            groupBy = groupBy,
            normMethod = normMethod,
            useGroups = useGroups,
            minCells = minCells,
            region = region,
            tileSize = tileSize,
            threads = threads,
            verbose = verbose,
            logFile = logFile
        )
    ArchR:::.logThis(split(df, df[, 3]), ".bulkTracks df", logFile = logFile)
    if (!is.null(ylim)) {
        ylim <- quantile(df$y, ylim)
        df$y[df$y < ylim[1]] <- ylim[1]
        df$y[df$y > ylim[2]] <- ylim[2] 
    }
    else {
        ylim <- c(0, max(df$y) * 1.03)

    }
    uniqueGroups <- gtools::mixedsort(unique(paste0(df$group)))
    if (!is.null(useGroups)) {
        uniqueGroups <- unique(useGroups)
    }
    df$group <- factor(df$group, levels = uniqueGroups)
    title <-
        paste0(as.character(seqnames(region)),
               ":",
               start(region) -
                   1,
               "-",
               end(region),
               " ",
               title)
    allGroups <-
        gtools::mixedsort(unique(
            getCellColData(
                ArchRProj = ArchRProj,
                select = groupBy,
                drop = TRUE
            )
        ))
    if (is.null(pal)) {
        pal <- suppressWarnings(paletteDiscrete(values = allGroups))
    }
    p <-
        ggplot(df, aes_string("x", "y", color = "group", fill = "group")) +
        geom_area(stat = "identity") + facet_wrap(facets = ~ group,
                                                  strip.position = "right",
                                                  ncol = 1) + ylab(
                                                      sprintf(
                                                          "Coverage\n(Norm. ATAC Signal Range (%s-%s) by %s)",
                                                          round(min(ylim), 2),
                                                          round(max(ylim), 2),
                                                          normMethod
                                                      )
                                                  ) +
        scale_color_manual(values = pal) + scale_fill_manual(values = pal) +
        scale_x_continuous(limits = c(start(region), end(region)),
                           expand = c(0, 0)) + scale_y_continuous(limits = ylim,
                                                                  expand = c(0, 0)) +
        theme_ArchR(
            baseSize = baseSize,
            baseRectSize = borderWidth,
            baseLineSize = tickWidth,
            legendPosition = "right",
            axisTickCm = 0.1
        ) + theme(
            panel.spacing = unit(0,
                                 "lines"),
            axis.title.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            strip.text = element_text(
                size = facetbaseSize,
                color = "black",
                margin = margin(0, 0.35, 0, 0.35,
                                "cm")
            ),
            strip.text.y = element_text(angle = 0),
            strip.background = element_rect(color = "black")
        ) +
        guides(fill = "none", colour = "none") + ggtitle(title)
    p
}
