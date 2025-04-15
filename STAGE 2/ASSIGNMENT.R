nhanes_data <- read.csv("https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/R/nhanes.csv")
head(nhanes_data)
names(nhanes_data)
str(nhanes_data)
summary(nhanes_data)
clean_data<-na.omit(nhanes_data)
clean_data$Age <- as.numeric(ifelse(clean_data$Age %in% c("", "NA"), NA, as.character(clean_data$Age)))

hist(clean_data$Age, main="Age distribution",col="yellow",xlab="Age",border="white")
hist(clean_data$BMI, main="BMI Distribution", col="skyblue")
hist(clean_data$Weight, main="Weight (kg) Distribution", col="pink")
hist(clean_data$Weight * 2.2, main="Weight (lb) Distribution", col="lightgreen")

mean(clean_data$Pulse, na.rm = TRUE)
min(clean_data$BPDia, na.rm = TRUE)
max(clean_data$BPDia, na.rm = TRUE)
library(ggplot2)

ggplot(clean_data, aes(x=Height, y=Weight, color=Gender)) +
  geom_point() +
  labs(title="Weight vs Height by Gender")
ggplot(clean_data, aes(x=Height, y=Weight, color=Diabetes)) +
  +     geom_point() +
  +     labs(title="Weight vs Height by Diabetes")
 
   ggplot(clean_data, aes(x=Height, y=Weight, color=SmokingStatus)) +
  +     geom_point() +
  +     labs(title="Weight vs Height by Smoking Status")

   # BMI by Diabetes
   t.test(BMI ~ Diabetes, data = nhanes_data)
   
   # Alcohol by Relationship Status
   t.test(AlcoholYear ~ RelationshipStatus, data =nhanes_data)
   
   # Age by Gender
    t.test(Age ~ Gender, data = nhanes_data)
   
