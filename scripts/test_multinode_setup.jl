using End2EndThermalImg
using Distributed

nw = End2EndThermalImg.setup_cluster_workers()

println("nworkers = ", nworkers())
println("workers = ", workers())
println("threads per worker = ", [remotecall_fetch(Threads.nthreads, p) for p in workers()])
println("hostname per worker = ", [remotecall_fetch(gethostname, p) for p in workers()])
