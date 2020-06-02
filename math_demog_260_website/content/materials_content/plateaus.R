library(tidyverse)

it <- read_table("./italy_bltcoh_1x1.txt")
it$Age <- as.numeric(it$Age)
it$mx <- as.numeric(it$mx)
it$qx <- as.numeric(it$qx)
it <- it %>% mutate(mx = if_else(mx == 0, as.numeric(NA), mx))
ggplot(it %>% filter(Age > 60, Year > 1900)) + 
  geom_line(aes(x = Age, y = log(mx), color = as.factor(Year))) + 
  ggtitle("log(mx) above age 60 after 1900")


# To detect plateaus above age 105 like the paper does:
# Run Gompertz on ages 45- 100 (there's some noise from the war) after 1900
# Gompertz: h(x) = a e^{b*x}
#         log(mx) = log(a) + bx
# We don't have hx so we'll have to use mx

gomp <- lm(log(mx) ~ Age, data = it %>% filter(Age > 45 & Age < 105 & Year > 1900))
new <- data.frame(Age = 105:109)
new$Gomp.pred <- exp(predict(gomp, new))
new <- new %>% left_join(it %>% filter(Age > 104, Year > 1900) %>% select(Age, Year,mx) %>% spread(Year, mx))
new <- new %>% gather("Year", "mx", -Age)
ggplot(new %>% filter(Year != "Gomp.pred")) + 
  geom_line(aes(Age, log(mx), color = Year)) + 
  geom_line(data = (new %>% filter(Year == "Gomp.pred")), aes(Age, log(mx) ), color = "black", size = 2) + 
  ggtitle("Gompertz prediction (black) ")

# Sort of!

# For men:
it_f <- read_table("./italy_mltcoh_1x1.txt")
it_f$Age <- as.numeric(it_f$Age)
it_f$mx <- as.numeric(it_f$mx)
it_f$qx <- as.numeric(it_f$qx)
it_f <- it_f %>% mutate(mx = if_else(mx == 0, as.numeric(NA), mx))
ggplot(it_f %>% filter(Age > 60 & Year > 1900)) + 
  geom_line(aes(x = Age, y = log(mx), color = as.factor(Year))) +
  ggtitle("Male log(mx) above age 60 after 1900 ") 

gomp <- lm(log(mx) ~ Age, data = it_f %>% filter(Age > 45 & Age < 95 & Year > 1900))
new <- data.frame(Age = 101:109)
new$Gomp.pred <- exp(predict(gomp, new))
new <- new %>% left_join(it_f %>% filter(Age > 95, Year > 1900) %>% select(Age, Year,mx) %>% spread(Year, mx))
new <- new %>% gather("Year", "mx", -Age)
ggplot(new %>% filter(Year != "Gomp.pred")) + 
  geom_line(aes(Age, log(mx), color = Year)) + 
  geom_line(data = new %>% filter(Year == "Gomp.pred"), aes(Age, log(mx) ), color = "black", size = 2) +
  ggtitle("Male: Gompertz Prediction (black)")


# Women
it_f <- read_table("./italy_fltcoh_1x1.txt")
it_f$Age <- as.numeric(it_f$Age)
it_f$mx <- as.numeric(it_f$mx)
it_f$qx <- as.numeric(it_f$qx)
it_f <- it_f %>% mutate(mx = if_else(mx == 0, as.numeric(NA), mx))
ggplot(it_f %>% filter(Age > 60 & Year > 1900)) + 
  geom_line(aes(x = Age, y = log(mx), color = as.factor(Year))) + 
  ggtitle("Female log(mx) above age 60 after 1900")

gomp <- lm(log(mx) ~ Age, data = it_f %>% filter(Age > 45 & Age < 95 & Year > 1900))
new <- data.frame(Age = 101:109)
new$Gomp.pred <- exp(predict(gomp, new))
new <- new %>% left_join(it_f %>% filter(Age > 95, Year > 1900) %>% select(Age, Year,mx) %>% spread(Year, mx))
new <- new %>% gather("Year", "mx", -Age)
ggplot(new %>% filter(Year != "Gomp.pred")) + 
  geom_line(aes(Age, log(mx), color = Year)) + 
  geom_line(data = new %>% filter(Year == "Gomp.pred"), aes(Age, log(mx) ), color = "black", size = 2) + 
  ggtitle("Female: Gompertz Prediction (Black) ")


