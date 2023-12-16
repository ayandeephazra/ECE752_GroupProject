The pursuit of efficient computer architecture has placed branch prediction at the forefront of processor performance optimization. Within this realm, the Tagged Geometric History Length (TAGE) branch predictor has been a notable benchmark since its introduction, especially for
its use of data compression principles derived from partial matching algorithms. This efficient approach has minimized pipeline stalls and maximized instruction throughput more effectively than previous models.

Despite TAGE’s prominence, its performance is capped by the static allocation of its table sizes, which fails to accommodate the diverse and dynamic nature of software workloads.

Our exploration into this domain underscores the need for a system that can adapt to varying workload demands, particularly those requiring longer history branches for accurate prediction. Our project endeavors to reproduce the dynamic branch prediction capabilities of TAGE, incorporating advancements from recent research. We are not innovating the reconfigurable architecture; instead, we are replicating the existing dynamic TAGE branch predictor within the gem5 simulation environment. Our objective is to verify the effectiveness of this established
method.
