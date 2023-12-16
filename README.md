The pursuit of efficient computer architecture has placed branch prediction at the forefront of processor performance optimization. Within this realm, the Tagged Geometric History Length (TAGE) branch predictor has been a notable benchmark since its introduction, especially for
its use of data compression principles derived from partial matching algorithms. This efficient approach has minimized pipeline stalls and maximized instruction throughput more effectively than previous models.

Despite TAGE’s prominence, its performance is capped by the static allocation of its table sizes, which fails to accommodate the diverse and dynamic nature of software workloads.

Our exploration into this domain underscores the need for a system that can adapt to varying workload demands, particularly those requiring longer history branches for accurate prediction. Our project endeavors to reproduce the dynamic branch prediction capabilities of TAGE, incorporating advancements from recent research. We are not innovating the reconfigurable architecture; instead, we are replicating the existing dynamic TAGE branch predictor within the gem5 simulation environment. Our objective is to verify the effectiveness of this established
method.

To streamline our project and make efficient progress, we decided to divide our efforts into two groups:
• Gem5 Team : This group primarily focused on software-based simulations using the Gem5 tool, such as -
– Gem5 Setup: Configure and set up Gem5 environments for simulating the Dynamic TAGE branch predictor

– Experiment Design: Design a comprehensive set of experiments to evaluate dynamic TAGE’s performance and accuracy under various SPEC workloads.

– Data Collection: Execute simulations, collect data, and analyze the results to assess dynamic TAGE’s behavior.

– Parameter Tuning: Fine-tune dynamic TAGE’s configuration parameters to optimize its prediction accuracy.

1

• Hardware Implementation Team: This group was responsible for the hardware-level implementation of the TAGE branch prediction algorithm in RTL.
– RTL Design: Develop the hardware architecture for dynamic TAGE
– Synthesis and Verification: Use hardware synthesis tools to convert RTL code into hardware descriptions and rigorously verify the correctness of the design.

– Testbench Creation: Create testbenches for RTL simulations to validate the hard-ware implementation’s functionality.

– Performance Evaluation: Analyze the hardware implementation in terms of resource utilization, power efficiency, and potential speedup in real processor architectures.
