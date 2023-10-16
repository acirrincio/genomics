#Function to return complementary DNA sequence from string 


getcomp <- function(x) {
    new <- ''
    for (i in 1:nchar(x)) {
      nuc <- substr(x,i,i)
      new <- str_c(new,ifelse(nuc=='A','T',
                              ifelse(nuc=='T','A',
                                     ifelse(nuc=='G','C',
                                            ifelse(nuc=='C','G','x')))))
      #print(new)
    }
      new
      return(new)
  }




# Example

seq <- 'ATGGCACAGCTTT'

getcomp(seq)


