#' Remove ms reps without required number of reps :/
#' 
#' @param data data.table
#' @param require.measurements minimum number of non 0 measurements
#' @param drop immediate drom measurements 
#' @param of set a flag
drop.zero.ms <- function(data, 
                         require.measurements = 2,
                         drop = FALSE,
                         flag = "drop_ms_rep")
{
    setkey(data, fragment_id, tech_id)
    data[,nzeros := sum(intensity > 0), by = list(fragment_id, tech_id)]
    data[, makenull := nzeros < require.measurements, by = list(fragment_id, tech_id)]
    data[,eval(flag) := FALSE]
    data[makenull == FALSE & intensity == 0, eval(flag) := TRUE]
    data[makenull == TRUE, intensity := 0]
    data[,c("makenull", "nzeros") := NULL]
    if (drop)
    {
        data = data[get(flag)==FALSE]
        data[, eval(flag) := NULL]
    }
    data
}

reconstruct.tech.single <- function(data)
{
    angDist <- function(data, f1, f2)
    {
        t1 <- data[intensity != 0 & fragment_id == f1, tech_id]
        t2 <- data[intensity != 0 & fragment_id == f2, tech_id]
        techs <- intersect(t1,t2)
        #CONSTANT
        if (length(techs) < 3)
            return(0)
        
        r <- angleDist(data[fragment_id == f1 & tech_id %in% techs, intensity],
                       data[fragment_id == f2 & tech_id %in% techs, intensity]) ## KEYBY

        r[!is.finite(r)] <- pi
        r <- (pi-abs(r))/pi
        r
    }
    
    setkey(data, tech_id, fragment_id)
    recover.tech <- function(fragment.id, tech.id, precursor.id, num)
    {
        #print(num)
        wt <- data[precursor_id == precursor.id, 
                  list(fragment_id, precursor_id, tech_id, intensity)]
        setkey(wt, fragment_id, tech_id)
        pivot_tech <- unique(wt[tech_id != tech.id & intensity != 0, tech_id])
        #pivot_int  <- wt[fragment_id == fragment.id & tech_id %in% pivot_tech, intensity]
        wt[, dist := angDist(wt, fragment.id, fragment_id), by=fragment_id]
        #wt[,dist:=1]
        pivot_mn <- data[fragment_id == fragment.id & tech_id %in% pivot_tech, mean(intensity)]
        
        wt[intensity >0, scale := pivot_mn/mean(intensity), by=fragment_id]
        i = wt[precursor_id     == eval(precursor.id) & 
               tech_id      == eval(tech.id) & 
               fragment_id  != eval(fragment.id) &
               intensity    != 0,
           c(sum(dist*intensity*scale)/sum(dist), .N)]        
        if (i[2] > 3)
            return(i[1])
        else
            return(0)
    }
    #browser()    

    data[, nzeros := sum(intensity > 0), by = fragment_id]
    #CONSTANT
    rec_candidates = unique(data[nzeros > 3 & intensity == 0,list(precursor_id),
                                 by=list(fragment_id, tech_id)])
    setkey(data, precursor_id, tech_id, fragment_id)
    
    rec_candidates[,r_intensity := recover.tech(fragment_id, tech_id, precursor_id, .I),
                   by = list(fragment_id, tech_id)]
    rec_candidates[, list(fragment_id, tech_id, r_intensity)]
}

reconstruct.tech.multiple <- function(data)
{
    angDist <- function(frame, f1, f2)
    {
        t1 <- frame[intensity != 0 & fragment_id == f1, tech_id]
        t2 <- frame[intensity != 0 & fragment_id == f2, tech_id]
        techs <- intersect(t1,t2)
        #CONSTANT
        if (length(techs) < 3)
            return(0)
        
        r <- angleDist(frame[fragment_id == f1 & tech_id %in% techs, intensity],
                       frame[fragment_id == f2 & tech_id %in% techs, intensity]) ## KEYBY
        
        r[!is.finite(r)] <- pi
        r <- (pi-abs(r))/pi
        r
    }
    
    setkey(data, tech_id, fragment_id)
    recover.tech <- function(fragment.id, tech.id, precursor.id, num)
    {
        #print(num)
        wt <- data[precursor_id == precursor.id, 
                  list(fragment_id, precursor_id, tech_id, intensity)]
        setkey(wt, fragment_id, tech_id)
        
        #wt[, dist := angDist(wt, fragment.id, fragment_id), by=fragment_id]
        
        wtech <- wt[tech_id == tech.id, list(fragment_id, ti = intensity)]
        setkey(wtech, fragment_id)
        setkey(wt, fragment_id)
        wt <- wt[wtech]
        wfrag <- wt[fragment_id == fragment.id, list(tech_id, fi = intensity)]
        setkey(wfrag, tech_id)
        setkey(wt, tech_id)
        wt <- wt[wfrag]
        wt[, est := ti*fi/intensity]
        
        wt <- wt[is.finite(est) & est >0]
        
        #Find mode
        if (length(wt$est) > 3)
        {
            dens <- density(wt$est)
            i <- dens$x[which.max(dens$y)]
        }else{
            i <- 0
        }
        
        #i = weighted.mean(wt$est, wt$dist)
        i

    }
    #browser()    
    
    data[, nzeros := sum(intensity > 0), by = fragment_id]
    #CONSTANT
    rec_candidates = unique(data[nzeros > 3 & intensity == 0,list(precursor_id),
                                 by=list(fragment_id, tech_id)])
    setkey(data, precursor_id, tech_id, fragment_id)
    
    rec_candidates[,r_intensity := recover.tech(fragment_id, tech_id, precursor_id, .I),
                   by=list(fragment_id, tech_id)]
    rec_candidates[,list(fragment_id, tech_id, r_intensity)]
}
