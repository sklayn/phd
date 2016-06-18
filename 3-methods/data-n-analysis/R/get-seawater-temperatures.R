library(XML)
# get the page source. It comes as a character vector contained in a single line:
url <- "http://www.stringmeteo.com/synop/sea_water.php?year=2012#sel"
raw.html <- readLines(url)

# set the encoding of the scraped page; otherwise Cyrillic characters (month and station names, date) get messed up during the parsing
Encoding(raw.html) <- "UTF-8"

# parse the tables (one month = one table)
scraped.tables <- htmlTreeParse(test, useInternalNodes = T, encoding = "UTF-8")
tables.list <- scraped.tables["//table"]

# use XPath expressions to extract and parse the values to obtain clean tables...

### GAVE UP AND COPY-PASTED THE DATA I NEEDED BY HAND; TO REDO AND FINISH THIS SOME FINE DAY WHEN MORE TIME / NERVES TO DEAL WITH HORRIFYING HTML TABLE AT SOURCE... 
