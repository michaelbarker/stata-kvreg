
/*******************************************************************************
 Set Moptimize Init Options Directly
 Primary and Secondary equation subcommands will call this to pass
 their equation-specific options.
*******************************************************************************/

program define np_moptimize_options
    version 11
    # delimit ;
    syntax name(name=M),
        [random repeat(integer 10)
         DIFficult TECHnique(string) 
         ITERate(numlist max=1 >0 integer) 
         TOLerance(numlist max=1 >0) LTOLerance(numlist max=1 >0) 
         NRTOLerance(numlist max=1 >0) NONRTOLerance 
         TRace GRADient showstep HESSian SHOWTOLerance noLOg 
         query initoptions(string asis)]
        ;
    # delimit cr
    
    if "`random'"=="random" {
        mata: moptimize_init_search_random(`M',"on")
        mata: moptimize_init_search_repeat(`M',`repeat')
    }

    if "`difficult'"=="difficult" {
        mata: moptimize_init_singularHmethod(`M',"hybrid")
    }

    if "`technique'"!="" {
        mata: moptimize_init_technique(`M', "`technique'")
    }

    if "`iterate'"!="" {
        mata:  moptimize_init_conv_maxiter(`M', `iterate')
    }
   
    if "`tolerance'"!="" {
        mata: moptimize_init_conv_ptol(`M', `tolerance') 
    }
    
    if "`ltolerance'"!="" {
        mata: moptimize_init_conv_vtol(`M', `ltolerance')
    }
    
    if "`nrtolerance'"!="" {
        mata: moptimize_init_conv_nrtol(`M', `nrtolerance')
    }

    if "`nonrtolerance'"=="nonrtolerance" {
        mata: moptimize_init_conv_ignorenrtol(`M', "on")
    }
   
    if "`trace'"=="trace" {
        mata: moptimize_init_trace_coefs(`M', "on")
    }
    
    if "`gradient'"=="gradient" {
        mata: moptimize_init_trace_gradient(`M', "on")
    }
   
    if "`showstep'"=="showstep" {
        mata: moptimize_init_trace_step(`M', "on")
    }
   
    if "`hessian'"=="hessian" {
        mata: moptimize_init_trace_Hessian(`M', "on")
    }
   
    if "`showtolerance'"=="showtolerance" {
        mata: moptimize_init_trace_tol(`M', "on")
    }
   
    if "`log'"=="nolog" {
        mata: moptimize_init_tracelevel(`M', "none")
    }

    gettoken option initoptions : initoptions , parse(", ") bind
    while `"`option'"' != "" {
        if `"`option'"' !="," {
            local option = subinstr(`"`option'"' , "(" , "(`M'," , 1)
            mata: moptimize_init_`option'
        }
        gettoken option initoptions : initoptions , parse(", ") bind
    }

    if "`query'"=="query" {
        mata: moptimize_query(`M')
    }

end


