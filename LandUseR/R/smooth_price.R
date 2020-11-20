

smooth_pr <- function() {
    x = fread(file.path(datadir(),"nico-output","FRA_model.csv"))
    z = x[!is.na(P_rural),.(year,P_rural)]
    plot(P_rural ~ year, z)
    lines(z$year,predict(lm(P_rural ~ ns(year,knots = c(1850,1950,1990)),z)))
}

