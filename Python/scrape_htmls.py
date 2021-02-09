#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 11:09:53 2021

@author: sohailnizam
"""

#-----------------------------------------------------------------------#
# Scrape cvR2, CI data from slapnap html and sample size from csv files #
#-----------------------------------------------------------------------#


#import necessary libraries
import os
from bs4 import BeautifulSoup
import pandas as pd
from glob import glob

#the bnabs which have directories in ic50ic80/docker_output 
#bnabs = ["vrc01", "10-1074", "10-996", "2f5", "2g12", "35o22",
#         "3bnc117", "4e10", "8anc195", "b12", "ch01", "hj16",
#         "nih45-46", "pg16", "pg9", "pgdm1400", "pgt121", 
#         "pgt128", "pgt135", "pgt145", "pgt151", "vrc-ch31",
#        "vrc-pg04", "vrc01", "vrc03", "vrc07", "vrc26.08",
#        "vrc26.25", "vrc29.03", "vrc34.01"]

bnabs = os.listdir("./ic50ic80/docker_output")


def scrape_tables(bnabs):
    
    '''
    This function iterates over the directories of the bnabs
    to extract info from both the .html report and the .csv data files.
    From the .html get the cross validated R2 as well as the CI lb and ubs
    for both ic50 and ic80. From the .csv file get the number of non NA values
    that bnab has for ic50 and for ic80. All extracted info stored in a pd df.
    
    Input: list of bnab name strings
    Output: pandas df with 2 rows per bnab (1 for ic50, 1 for ic80)
            and 6 columns (bnab, method, r2, cil, ciu, n)
    '''
    #init empty list to hold scraped r2s, lbs, and ubs
    cvR2_data = []

    
    for bnab in bnabs:
        
        #assume ic50ic80 is the current working dir
        path = "./docker_output/"
        path += bnab
        #report_path = path + "/report_" + bnab.upper() + "_" + date + ".html"
        
        
        #if the container failed, the dir does not have a report or csv file
        #also skips the .gitkeep which may mistakenly be put in the nab list
        try:
            #pick out the file that starts with 'report' ends with '.html'
            #this method is not date specific
            report_path = glob(path + "/report*.html")[0]
            report = open(report_path).read()
            
            #pick out the file that starts with 'slapnap' ends with '.csv'
            #this method is not date specific
            data_path = glob(path + "/slapnap*.csv")[0]
            bnab_data = pd.read_csv(data_path)
            
        except IndexError:
            continue
        
        #scrape the html with beautiful soup
        soup = BeautifulSoup(report)
        #extract the second table in the html
        cvR2_table = soup.findAll("table")[1]
        
        #init a list to hold the data for this bnab
        output_rows = []
        #two rows not including header: ic50 and ic80
        for row_num, table_row in enumerate(cvR2_table.findAll("tr")[1:]):
            
            #three columns not including label
            columns = table_row.findAll('td')
            #init one row of the data
            output_row = []
            
            #first append the method (ic50 or ic80)
            if row_num == 0:
                output_row.append('ic50')
            else:
                output_row.append('ic80')
            
            
            #next append the bnab name
            output_row.append(bnab)
            
            #finally append the info from the 3 columns of this row 
            #(r2, cil, ciu)
            for column in columns[1:]:
            
                output_row.append(column.text)
            
            #in total, output_rows will be a list of two lists
            #each sublist with entries: method, bnab, r2, cil, ciu
            #first sublist for ic50, second sublist for ic80
            output_rows.append(output_row)
            
        #now extract the sample sizes for this bnab from the csv
        n_ic50 = bnab_data['ic50'].count()
        n_ic80 = bnab_data['ic80'].count()
        
        #append the sample sizes to the appropriate sublist
        output_rows[0].append(n_ic50)
        output_rows[1].append(n_ic80)
        
        #add the rows for this bnab to the running list
        cvR2_data += output_rows
        #print(bnab + " done")
    
    
    #cast the data list of lists to pandas df
    cvR2_dataframe = pd.DataFrame(cvR2_data, columns = ["method", "bnab", "r2", 
                                                        "cil", "ciu", "n"])

    #return sorted to ic50 and ic80 are separate
    return(cvR2_dataframe.sort_values(by = ['method', 'bnab']))



#call the function on the bnabs
result_df = scrape_tables(bnabs)
        
#write to csv
result_df.to_csv("./Python/cvr2_df.csv")
 

    
    
    