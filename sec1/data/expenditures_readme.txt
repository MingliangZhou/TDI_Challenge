House Office Expenditure Data
Last updated November 28, 2016.
https://projects.propublica.org/represent/expenditures

Members of the House of Representatives get an annual budget for their Washington and district offices, but how they spend it is up to them. There are some rules: It can’t be used for personal or campaign expenses, and there is no reserve source of money if lawmakers spend all of their allowances.

Lawmakers also are required to report the recipients of their office spending, and since 2009 the Sunlight Foundation has been taking the PDF files published by the House and converting them into text files useful for analysis and research. As of November 2016, ProPublica has taken over both the collection and hosting of these files. They can be examined using spreadsheet or database software. 

How We Collect This Data

Each quarter we take the report published by the House and generate two text files: One contains summary information for each office and category of spending (some examples include “Personnel Compensation” and “Travel”), and the other contains details of each recipient of office spending and its purpose. Note that the data has not been standardized (meaning that "AT&T" might also appear as "A.T.&T."), so simple aggregation on the recipient could result in multiple totals for the same individual or entity, depending on the spelling. Individual recipients can be paid by more than one office or lawmaker in some cases.

Most of the records are connected to lawmaker offices, but the files also contain spending records for House committees and administrative offices, in addition to leadership organizations such as the Speaker of the House and the two parties' leaders.

Before you dig into the data to find out how the House spends its money, check out this slideshow from a Sunlight training webinar on how to navigate the database. You should also read this post that explains discrepancies with how the House reports lawmakers' spending, and gives guidelines on how to use the data. For some examples of how this data can be used, check out Sunlight's analyses of congressional capacity during reduced operating budgets and staff turnover among members' offices.


Data Dictionary

Summary files

BIOGUIDE_ID – the official ID of members of the House (http://bioguide.congress.gov/biosearch/biosearch.asp)
OFFICE – the name of the House office
YEAR – the calendar year
QUARTER – the quarter of the year
CATEGORY – broad description of spending
YTD – year to date amount spent by office in that category 
AMOUNT – amount spent by office in that category in quarter

Detail files

Has BIOGUIDE, OFFICE, QUARTER, YEAR, CATEGORY, AMOUNT, plus:

PAYEE – name of recipient
PURPOSE – specific purpose of spending
DATE -  date of payment (optional)
START DATE – beginning of period which payment covers
END DATE – end of period which payment covers
TRANSCODE – House transaction code
TRANSCODELONG – description of House transaction code
RECORDID – House record number
RECIP (orig.) - original (non standardized) recipient
