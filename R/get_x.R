get_x = function(formula, meta_data) {
    opt = options(na.action = "na.pass") # Keep NA's in rows of x
    on.exit(options(opt)) # Switch it back
    x = model.matrix(formula(paste0("~", formula)), data = meta_data)
    return(x)
}
