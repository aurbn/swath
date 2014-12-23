reconstruct.measurements <- function(data, completeness = 5)
{
   # browser()
    # Find all "zero" measurements
    data[,zero := intensity <= 0]
    data[,zeros := sum(zero), by = fragment_id]
    data[,rec_target := zeros != 0 & zeros < completeness]
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
        rec_runs = wt[fragment_id == f & zero == TRUE, run_id]
        for (r in rec_runs)
        {
            pivot_runs <- wt[fragment_id == f & zero == FALSE, run_id]
            if (length(pivot_runs) < 2)
                next
            pivot_int  <- wt[fragment_id == f & run_id %in% pivot_runs,intensity]
            
            #Angular distance
            ad <- make_ad(pivot_int)
            dist <- wt[run_id %in% pivot_runs, ad(intensity), by=fragment_id]$V1
            
            Is = c()
            for (pr in pivot_runs)
            {
                I = pivot_int*wt[run_id != get(pr), ]
            }
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
                ddist <- rep(1, length(ddist))
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

reconstruct.measurements.simple <- function(data, completeness = 5, recreq = 3)
{
    # browser()
    # Find all "zero" measurements
    data[,zero := intensity <= 0]
    data[,zeros := sum(zero), by = fragment_id]
    data[,rec_target := zeros != 0 & zeros < completeness]
    data[,recovered := FALSE]
    
    make_ad <- function(pivot)
    {
        function(int)
        {
            good = pivot!=0 & int != 0 
            #CONSTANT
            if (sum(good) < 5)
                return(0)
            r <- angleDist(pivot[good], int[good])
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
        
        # "Current frane"
        setkey(data, precursor_id)
        wt = data[precursor_id == prec]
        pivot_int  <- wt[fragment_id == f,intensity]
        #Angular distance
        
        ad <- make_ad(pivot_int)
        wt[, dist := ad(intensity), by=fragment_id]
        
        #Between target and other fragment
        #Scaling fator (for different absolute intensisies)
        wt[intensity >0, scale := mean(pivot_int[pivot_int>0])/mean(intensity), by=fragment_id]
        
        #Recovery _target_measurements_
        rec_runs = wt[fragment_id == f & zero == TRUE, run_id]
        for (r in rec_runs)
        {
            i <- wt[precursor_id == eval(prec) & 
                    run_id       == eval(r) & 
                    fragment_id  != eval(f) &
                    zero         != TRUE,
                    sum(dist*intensity*scale)/sum(dist)]
            
            if (nrow(rest) < recreq){
                next
            }else{
                data[run_id == r & fragment_id == f, intensity := i]
                data[run_id == r & fragment_id == f, recovered := TRUE]
            }
        }
    }
    return(data)
}