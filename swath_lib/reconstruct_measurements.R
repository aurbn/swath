reconstruct.measurements <- function(data, completenbess = 20)
{
   # browser()
    # Find all "zero" measurements
    data[,zero := intensity <= 0]
    data[,zeros := sum(zero), by = fragment_id]
    data[,rec_target := zeros != 0 & zeros < 10]
    data[,recovered := FALSE]

    make_ad <- function(pivot)
    {
        function(int)
        {
            r <- angleDist(pivot, int)
            r[!is.finite(r)] <- pi
            r <- (pi-abs(r))/pi
            r
        }
    }
    
    fragments <- unique(data[rec_target == TRUE, fragment_id])
    n.fragments = length(fragments)
    i.fragments = 1
    
    create.directory(subDir = "refs")
    
    for (f in fragments)
    {
        loginfo("Recontruction %i of %i fragment, %.1f%%",
                i.fragments, n.fragments, i.fragments*100/n.fragments)
        i.fragments = i.fragments + 1
        
        prec = unique(data[fragment_id==f, precursor_id])
        if (length(prec)!=1)
        {
            stop("More than one precursor for fragment!")
        }
        wt = data[precursor_id == prec]
        
        #Recovery _target_measurements_
        recs = wt[fragment_id == f & zero == TRUE, run_id]
        for (r in recs)
        {
            pivot_runs <- wt[fragment_id == f & zero == FALSE, run_id]
            if (length(pivot_runs) < 2)
                next
            pivot_int  <- wt[fragment_id == f & run_id %in% pivot_runs,intensity]
            
            #Angular distance
            ad <- make_ad(pivot_int)
            dist <- wt[run_id %in% pivot_runs, ad(intensity), by=fragment_id]$V1
            
            Is = c()    
            for (rr in pivot_runs)
            {
                I1 = wt[run_id == r, intensity]
                I2 = wt[run_id == rr,intensity]
                Iok = I2 > 0
                if (sum(Iok) < 2)
                    next
                ddist = dist[Iok]
                I1 = I1[Iok]
                I2 = I2[Iok]
                I = wt[fragment_id == f & run_id == rr, intensity]*
                    sum(na.omit(I1*ddist/I2))/sum(dist)
                Is = c(Is, I)
            }
            png(paste(sep="", "refs/", i.fragments, "_", r))
            plot(density(Is), main = paste(sep="_", mean(Is), length(Is)))
            dev.off()
            Is = mean(Is)
            data[run_id == r & fragment_id == f, intensity := Is]
            data[run_id == r & fragment_id == f, recovered := TRUE]
        }
    }

}