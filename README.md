# R-code-for-epimedic-models
This profile includes all the relative data:
    Y: weekly reported Measles epidemics data in Liverpool, Birmingham and London between year 1949-1953
    B: the real birth data for year 1945-1949 and 
    N: population data for year 1949-1953 for each city 
    D: the distances between each pair of cities
We try to use both standard SMC method and GIRF method introduced by Joonha Park and Edward (2017) to simulate particle filters and do paramter inference using R.
Due to the expensive computaional cost, the number of particel filters J used here are rather small but we can explore the performance of these two methods as well.
