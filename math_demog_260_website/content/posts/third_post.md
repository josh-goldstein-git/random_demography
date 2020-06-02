---
title: "Increasing inequality?"
date: 2020-02-13T20:15:51-08:00
draft: true
---

The Waldron [analysis](https://www.ssa.gov/policy/docs/ssb/v67n3/v67n3p1.html) presents the following table. 

{{< figure src="/posts/waldron.png" title="">}}

If we simulate two heterogeneous groups, one which always has twice
the baseline mortality level as the other, then we will see the
observed hazard ratio converge with age as mortality selectively
removes the most frail individuals. If we let the baseline levels
of both groups decline steadily, simulating the process of mortality
improvement, then this will make it appar that the hazard ratio is
rising over time.

Here is the results of such a simulation:

| cohort | age 60  | age 70 | age 80 | age 90 | age 100|
| ------ | ---:|---:|---:| ---:| ---:|
| 1      |1.44 |1.23|1.10| 1.04|1.01|
| 2      |1.51 |1.27|1.12| 1.05|  NA|
| 3      |1.57 |1.33|1.15|   NA|  NA|
| 4      |1.63 |1.38|  NA|   NA|  NA|
| 5      |1.68 |  NA|  NA|   NA|  NA|

The table shows the same pattern as Waldron's estimates, however the
magnitude of change is considerably smaller. 

One could try calibrating the simulation to match observed mortality
rates in the United States for the time period Waldron considers.

The code used to generate the simulation is [here](../../materials_content/waldron_simu.R)
