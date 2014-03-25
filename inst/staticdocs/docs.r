library(staticdocs)
library(grid)
list(
    readme = "",
    
    index = list(
        sd_section("Analysis functions",
            "These functions are used for analysis",
            c(
                'dGSEA',
                'dGSEAview',
                'dGSEAwrite',
                'dPvalAggregate',
                'dNetInduce',
                'dBUMfit',
                'dBUMscore',
                'dSVDsignif',
                'dFDRscore',
                'dNetFind',
                'dNetPipeline',
                'dNetReorder',
                'dCommSignif',
                'dNetConfidence',
                'dContrast',
                'dRWR'
            )
        ),
        sd_section("Visualisation functions",
            "These functions are used for visualisation",
            c(
                'visGSEA',
                'visNet',
                'visNetMul',
                'visNetReorder',
                'visNetArc',
                'visNetCircle',
                'visColormap',
                'visColoralpha',
                'visHeatmap',
                'visHeatmapAdv',
                'visTreeBootstrap'
            )
        ),
        sd_section("Built-in databases in human",
            "",
            c(
                "org.Hs.eg",
                "org.Hs.egDO", 
                "org.Hs.egGOBP",
                "org.Hs.egGOCC",
                "org.Hs.egGOMF",
                "org.Hs.egHPMI", 
                "org.Hs.egHPON",
                "org.Hs.egHPPA",
                "org.Hs.egMP",
                "org.Hs.egPS",
                "org.Hs.string"
            )
        ),
        sd_section("Built-in databases in mouse",
            "",
            c(
                "org.Mm.eg",
                "org.Mm.egDO", 
                "org.Mm.egGOBP",
                "org.Mm.egGOCC",
                "org.Mm.egGOMF",
                "org.Mm.egHPMI", 
                "org.Mm.egHPON",
                "org.Mm.egHPPA",
                "org.Mm.egMP",
                "org.Mm.egPS",
                "org.Mm.string"
            )
        ),
        sd_section("Built-in databases in arabidopsis",
            "",
            c(
                "org.At.string"
            )
        ),
        sd_section("Built-in datasets",
            "",
            c(
                "Hiratani_TableS1",
                "CLL"
            )
        )
    ),
    
    if(0){
    icons = list(  
        eCal = sd_icon({
          textGrob("Common", rot = 45, gp = gpar(cex = 1))
        }),
        visRunES = sd_icon({
          textGrob("Hot", rot = 45, gp = gpar(cex = 1.2))
        }),
        eView = sd_icon(inherit = "eCal"),
        visNet = sd_icon(inherit = "visRunES")
    )
    }
)
