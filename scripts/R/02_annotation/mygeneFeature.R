mygeneFeature = function(bam, type = "within"){
              single.bam.gr = get(bam)
              single.bam.gr$region = NULL
              single.bam.new.gr = test.bam.gr
              single.bam.new.gr = test.bam.new.gr[-(1:length(test.bam.new.gr))]
              single.other.gr = NULL

              ##sense
              ## type = "within": reads fully contained in the feature (default).
              ## type = "union": either the read is fully contained in the feature
              ## OR the feature is fully contained in the read (union containment).
              ## type = "any": any overlap (>= 1 bp) between read and feature.
              for(m in 1:length(geneFeature.list)){
                    geneFeature.name = geneFeature.list[m]
                    geneFeature.id = gsub("refGene\\.NM\\.", "", geneFeature.name)
                    geneFeature.id = gsub("\\.grl", "", geneFeature.id)
		    geneFeature.id = gsub("\\.gr", "", geneFeature.id)

                    single.bam.ol = NULL
                    single.bam.hit.gr = NULL
                    if(type == "union"){
                      ol1 = findOverlaps(single.bam.gr, get(geneFeature.name), type = "within")
                      ol2 = findOverlaps(get(geneFeature.name), single.bam.gr, type = "within")
                      hit.ids = unique(c(queryHits(ol1), subjectHits(ol2)))
                    }else if(type == "any"){
                      single.bam.ol = findOverlaps(single.bam.gr, get(geneFeature.name), type = "any")
                      hit.ids = unique(queryHits(single.bam.ol))
                    }else{
                      single.bam.ol = findOverlaps(single.bam.gr, get(geneFeature.name), type = "within")
                      hit.ids = unique(queryHits(single.bam.ol))
                    }
                    single.bam.hit.gr = single.bam.gr[hit.ids]
                    single.bam.hit.gr$region = rep(geneFeature.id, length(single.bam.hit.gr))
                    single.bam.new.gr = c(single.bam.new.gr, single.bam.hit.gr)
                    if(length(hit.ids) == 0){
                     single.bam.gr = single.bam.gr
                    }else{
                     single.bam.gr = single.bam.gr[!single.bam.gr$rm.key %in% single.bam.gr$rm.key[hit.ids]]
                    }
              }
             for(g in 1:length(geneFeature.list)){
                    geneFeature.name = geneFeature.list[g]
                    geneFeature.id = gsub("refGene\\.NM\\.", "", geneFeature.name)
                    geneFeature.id = gsub("\\.grl", "", geneFeature.id)
		    geneFeature.id = gsub("\\.gr", "", geneFeature.id)

                    single.bam.ol = NULL
                    single.bam.hit.gr = NULL
                    single.bam.ol = findOverlaps(single.bam.gr, get(geneFeature.name), type = "any", ignore.strand = T)
                    single.bam.hit.gr = single.bam.gr[unique(queryHits(single.bam.ol))]
                    single.bam.hit.gr$region = rep(paste("AS.", geneFeature.id, sep = ""), length(single.bam.hit.gr))
                    single.bam.new.gr = c(single.bam.new.gr, single.bam.hit.gr)
                    if(length(unique(queryHits(single.bam.ol))) == 0){
                     single.bam.gr = single.bam.gr
                    }else{
                     single.bam.gr = single.bam.gr[!single.bam.gr$rm.key %in% single.bam.gr$rm.key[unique(queryHits(single.bam.ol))]]
                    }
              }

              single.bam.other.gr = single.bam.gr
              single.bam.other.gr$region = rep("intergenic", length(single.bam.other.gr))
              single.bam.new.gr = c(single.bam.new.gr, single.bam.other.gr)
              return(single.bam.new.gr)
}
