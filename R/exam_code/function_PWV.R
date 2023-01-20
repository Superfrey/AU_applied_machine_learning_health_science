### Coleset two means


closest_mean <- function(v) {
    val_count <- 3 - (is.na(v[1]) + is.na(v[2]) + is.na(v[3]))

    if (val_count==0) {
        x <- NA
    } else if (val_count %in% c(1,2)) {
        x <- mean(v, na.rm = TRUE)
    } else if (val_count==3) {
        if (abs(v[2]-v[3]) <= abs(v[1]-v[3]) & abs(v[2]-v[3]) <= abs(v[1]-v[2])) x <- mean(c(v[2],v[3]))
        if (abs(v[1]-v[3]) <= abs(v[1]-v[2]) & abs(v[1]-v[3]) <= abs(v[2]-v[3])) x <- mean(c(v[1],v[3]))
        if (abs(v[1]-v[2]) <= abs(v[1]-v[3]) & abs(v[1]-v[2]) <= abs(v[2]-v[3])) x <- mean(c(v[1],v[2]))
    }

    return(x)
}
