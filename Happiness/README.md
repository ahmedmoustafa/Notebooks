# World Happiness Report

![Smile](http://inspiration.rehlat.com/wp-content/uploads/2018/03/World-Happiness-Report-2018.jpg)

[The World Happiness Report](https://worldhappiness.report/) is a landmark survey of the state of global happiness that ranks countries by how happy their citizens perceive themselves to be.

The 2019 dataset includes the following columns:

| Column | Description |
| -- | -- |
| `country` | Country or Region |
| `category` | UN classification: Developed, Developing, Underdeveloped, and Transiting |
| `score` | A metric measured by asking the sampled people the question: "How would you rate your happiness on a scale of 0 to 10 where 10 is the happiest." |
| `gdp_per_capita` | GDP per capita |
| `social_support` | Social support |
| `healthy_life_expectancy` | Healthy life expectancy |
| `freedom_to_make_life_choices` | Freedom to make life choices |
| `generosity` | Generosity |
| `perceptions_of_corruption` | Perceptions of corruption |


Using R analysis, answer the following questions:

1. What is the happiness `score` for `Egypt`? **Hint**: use `R` code to find it :wink:

2. How many countries do have a happiness `score` > `Egypt`'s `score`?

3. How many countries do have a happiness `score` < `Egypt`'s `score`?

4. What are the medians of the `score` for each `category`? Compute the medians, then plot the values.

5. Is the happiness `score` correlated with `gdp_per_capita`? Show visually and numerically. **Hint**: for the numerical part, use the [`cor`](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/cor) and [`cor.test`](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/cor.test) functions

6. Other than the `Developed` countries, which countries do have happiness `score` > median of the scores + standard deviation of the scores? **Hint**: the standard deviation can be calculated using the [`sd`](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/sd) function

7. Write `R` code to generate the figure below. How would you interpret the figure?
![Generosity](generosity.png)

8. Does the happiness `score` change significantly with the country `category`? Show visually and numerically. **Hint**: for the numerical part, start with the [`aov`](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/aov) function

9. Which of the reported parameters has the most significant association with the `score`? **Hint**: exclude `country` and `category` from this part.

10. Which of the reported parameters has the largest effect on the `score`? **Hint**: exclude `country` and `category` from this part.

To start answering, copy [this notebook](Happiness.ipynb) into your Colaboratory on Google Drive.

Good luck :star2:
