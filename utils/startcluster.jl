hostname = readlines(`hostname`)[1] 
cmd = readlines(pipeline(`cat /proc/$(getpid())/cmdline`,`xargs -0 echo`))[1]

if hostname=="ux305fa"
    ### Sync the code to hpc for example with rsync; also recreate empty folders; removes unneded files as well
    run(`rsync -a --delete-excluded '--include=*.jl' '--include=*/' '--exclude=*' . 'hpc05:~/data'`)
    ## ssh -c to run the same command we used
    run(`ssh hpc05 -t "export PATH=/home/jerdmanis/julia/bin:\$PATH; cd ~/data; $cmd"`)
    ### sync back created jld files
    run(`rsync -a '--include=*.jld' '--include=*/' '--exclude=*' 'hpc05:~/data/' .`)
    info("Sync finished")
    exit()
elseif hostname=="hpc05"
    ### remove all jld files
    using ClusterManagers

    run(`mkdir -p data`)
    run(`date`)
    println("How many cores you will need?")
    np = parse(Int,readline())
    println("How much time will you need hh:mm:ss?")
    walltime = readline()
    if walltime==""
        walltime="01:00:00"
    end
    println("How much memory will you need GB?")
    memory = readline()
    if memory==""
        memory="2"
    end
    
    #addprocs_pbs(np,queue="q1",res_list="walltime=$walltime;mem=$(memory)GB")
    addprocs_pbs(np,queue="q1",res_list="walltime=$walltime,mem=$(memory)GB")

    ### here I also could redefine println and pmap for a better laod balancing
    @everywhere import Base.println
    @everywhere println(x::AbstractString) = run(`echo "$x"`)

    @everywhere BLAS.set_num_threads(1)
    
    info("Cluster Initialized")
    run(`date`)
end
