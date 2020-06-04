# Searched 31/05/2020
# terms: "covid"[Title] OR "sars-cov-2"[Title] OR "2019-ncov"[Title] OR "novel coronavirus"[Title] 
# 16,003 results

#install.packages(c("tidyverse", "lubridate", "gridExtra", "cowplot", "ggrepel",
#                 "ggpubr", "scales", "tidytext", "wordcloud"))

library("tidyverse")
library("lubridate")
library("gridExtra")
library("cowplot")
library("ggrepel")
library("ggpubr")
library("scales")
library("tidytext")
library("wordcloud")

# Load and prepare data -------------------------------------------------------

# Parse csv file from PubMed
rawdf <- as_tibble(read.csv(file = "covid_pubmed_200531_merged.csv", stringsAsFactors = FALSE)) %>%
  mutate_all(~ na_if(., "")) #replace blank cells with NA
rawdf

# Check dates
rawdf <- rawdf %>% mutate(Create.Date = as.Date(Create.Date, "%d/%m/%Y"))
sum(is.na(rawdf$Create.Date)) #0
summary(rawdf$Create.Date) #min = 2003-04-12
earlydf <- rawdf %>% filter(Create.Date < as.Date("2020-01-01"))
##View(earlydf) #these are related to other 'novel coronaviruses', so limit to 2020-01-01
rawdf <- rawdf %>% filter(Create.Date >= as.Date("2020-01-01"))
summary(rawdf$Create.Date) #first publication on 2020-01-19

# Remove duplicated titles and unnessesary columns
df <- rawdf %>%
  distinct(Title, .keep_all = TRUE) %>%
  select(Title, Create.Date, Journal.Book) %>%
  rename(title = "Title", date = "Create.Date", journal = "Journal.Book")
#n.b. .keep_all = T refers to keeping all columns, not rows


# Dates and counts ------------------------------------------------------------

# create sequence of dates for each day
dateseq <- seq.Date(from = as.Date("2020-01-01"), to = max(df$date), by = 1)
# create simplified tbl
datedf <- data.frame(date = dateseq, n = NA)
# for loop to get counts
for (i in 1:nrow(datedf)){
  datedf[i,2] <- sum(df$date == datedf[i,1])
}
datedf <- as_tibble(datedf)
datedf
# add cumulative counts
datedf <- datedf %>% mutate(cumn = cumsum(n))
datedf

# Plot of publications per day ------------------------------------------------

# helper function to split month into two (for x-axis labels)
bimonthly <- function(x) {
  x_range <- range(x, na.rm = TRUE)
  
  date_range <- c(
    floor_date(x_range[1], "month"),
    ceiling_date(x_range[2], "month")
  )
  monthly <- seq(date_range[1], date_range[2], by = "1 month")
  
  sort(c(monthly, monthly + days(14)))
}

#plot
plot_day <- datedf %>%
  ggplot(aes(x = date, y = n)) +
  geom_point(colour = "#0072B2") +
  stat_smooth(method = "loess", se = FALSE, colour = "black", linetype = "solid") +
  scale_y_continuous(trans = "log10") +
  scale_x_date(date_labels = "%d %b %y", breaks = bimonthly) +
  labs(x = "Date", y = "Number of publications per day") +
  theme_pubr(x.text.angle = 45) +
  NULL
plot_day

# Plot of cumulative publications
plot_cumulative <- datedf %>%
  ggplot(aes(x = date, y = cumn)) +
  geom_point(colour = "#0072B2") +
  geom_smooth(method = "loess", se = FALSE, colour = "black", linetype = "solid") +
  scale_y_continuous(trans = "log10") +
  scale_x_date(date_labels = "%d %b %y", breaks = bimonthly) +
  labs(x = "Date", y = "Cumulative number of publications") +
  theme_pubr(x.text.angle = 45) +
  NULL
plot_cumulative

pubs_grid <- plot_grid(plot_day, plot_cumulative, nrow = 1, labels = c("A", "B"))
pubs_grid
ggsave("figure1.png", pubs_grid, width = 14, height = 8)


#Summary statistics
#median publications per day
datedf %>% filter(n > 0) %>% summary() #since first publication


# Text mining -----------------------------------------------------------------

# create a column for month of publication
tidy_titles <- df %>%
  mutate(month = factor(case_when(
           grepl("-01-", date) ~ "January",
           grepl("-02-", date) ~ "February",
           grepl("-03-", date) ~ "March",
           grepl("-04-", date) ~ "April",
           grepl("-05-", date) ~ "May")),
        month = fct_relevel(month, "January", "February", "March", "April", "May")
  )

# keep track of publication number
tidy_titles <- tidy_titles %>%
  group_by(month) %>%
  mutate(pubnumber = row_number()) %>%
  ungroup()

# tidy data frame of title words
month_words <- tidy_titles %>%
  mutate(title = str_replace_all(title, "-", "_")) %>% #prevent hyphenated words being separated
  unnest_tokens(word, title) %>%
  mutate(word = str_replace_all(word, "_", "-")) %>% #put hyphens back
  count(month, word, sort = TRUE)

# summary of total title words each month
total_words <- month_words %>%
  group_by(month) %>%
  summarise(total = sum(n))

# join these
month_words <- left_join(month_words, total_words)

# some uninformative words to ignore
mystopwords <- tibble(word = c("2019-new", "covid-19", "sars-cov-2", "covid", "44", "000", "is"))

# remove mystopwords and add tf-idf metric for each word
month_words <- month_words %>%
  anti_join(mystopwords, by = "word") %>%
  bind_tf_idf(word, month, n) %>%
  arrange(desc(tf_idf))

month_words

# plot top 15 tf-idf words for each month
tf_idf_plot <-
  month_words %>%
  group_by(month) %>%
  top_n(15, tf_idf) %>%
  ungroup () %>%
  mutate(word = reorder_within(word, tf_idf, month)) %>%
  ggplot(aes(word, tf_idf, fill = month)) +
  geom_col(show.legend = FALSE) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = NULL, y = "tf-idf") +
  facet_wrap(~month, nrow = 2, scales = "free") +
  coord_flip() +
  scale_x_reordered() +
  theme_pubr() +
  NULL

ggsave("figure3.png", tf_idf_plot, width = 12, height = 6)


# Word cloud ------------------------------------------------------------------

# load stop words data from tidytext
data(stop_words)
mystopwords2 <- tibble(word = c("coronavirus", "2019-ncov", "covid-19", "sars-cov-2", "sars", "covid"))

# create tidy data frame for word cloud
word_cloud_df <- tidy_titles %>%
  mutate(title = str_replace_all(title, "-", "_")) %>% #prevent hyphenated words being separated
  unnest_tokens(word, title) %>%
  mutate(word = str_replace_all(word, "_", "-")) %>% #put hyphens back
  filter(!grepl("^[0-9.-]*$", word)) %>% # remove 'words' that are just digits +/- hyphen/dot
  anti_join(mystopwords2, by = "word") %>%
  anti_join(stop_words, by = "word") %>%
  count(word, sort = TRUE)
  
# plot word cloud (100 words max)
pal <- brewer.pal(8,"Dark2")
png("figure2.png", width=12, height=8, units='in', res=300)
set.seed(123)
wordcloud(word_cloud_df$word, word_cloud_df$n, random.order = FALSE, max.words = 100, colors=pal)
dev.off()
