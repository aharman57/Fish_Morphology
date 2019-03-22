# Fish_Morphology
Final Project for Biology 708

Data:
Our dataset consists of information on several morphological measurements of Lake Whitefish hatchlings that were reared in different water temperatures. The Age represents the number of days after hatching at which morpholigcal measurements were taken. Treatment values represent the temperature (in degrees Celsius) of the water in which eggs were incubated. After the icubation treatment (post-hatch) larvae were all kept at the same temperature (8 degC). Spawning adults were caught in Lake Huron and stripped of eggs/milt, and fertilized embryos brought back to the lab. Sampling was lethal, so each row represents a different fish (different fish were sampled at different ages). This experiment was repeated several years in a row. This is only 1 year of the dataset, so if my project is chosen I will acquire several more years of data. All the response variables are correlated, as they all come from the same individual. As the fish grows, so do all the other body features that we measure (except the yolk, which diminishes). A simple linear model won't account for the fact that these variables are co-related, so I want to be able to fit a model that accounts for the correlations. 

I want to determine the relative importance of each response, hierarchy of influence on length. Which variable is the best indicator of total length? I also want to be able to accurately predict total length based on a few morphological measurements, to make data collection easier.


Scientific Questions:
1. What is the influence of temperature on length and biomass? How does temperature influence the relationship between length and biomass?

2. What is the relationship between length and other morphological characteristics, and how does temperature affect these relationships?

3. What is the effect of temperature on 'yolk efficiency' (i.e. relationship between increase in length or biomass and reduction in yolk size over time)?

Analysis Plans:
Correlation tests to see how characters are correlated
Regression to test effects - temperature is continuous but not quite because interval scale (can't multiply/divide these numbers in a meaningful way)

Questions/Issues:
1. How do we deal with correlation between different morphological characteristics? Correlation makes it difficult to tell whether temperature is directly affecting a certain morphological feature, or if this is just due to correlation with another morphological feature.

2. How do we test interactions between multiple morphological characteristics so that we can understand if temperature influences the relationships between these characters?

3. Can we create a model that describes the relationships between different morphological characteristics so that we can infer the response of one character based on the data about other characters? We want to be able to accurately predict total length based on a few morphological measurements, to make data collection easier.