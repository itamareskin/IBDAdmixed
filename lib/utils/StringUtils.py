'''
Created on Dec 11, 2011

@author: itamares
'''

from libc.stdlib cimport malloc, free 
from libc.string cimport strlen, strcpy
import re

def multiple_replace(d, text): 
        """ Replace in 'text' all occurences of any key in the given
        dictionary by its corresponding value.    Returns the new tring.""" 
    
        # Create a regular expression    from the dictionary keys
        regex = re.compile("(%s)" % "|".join(map(re.escape, d.keys())))
    
        # For each match, look-up corresponding value in dictionary
        return regex.sub(lambda mo: d[mo.string[mo.start():mo.end()]], text)
    
class StringUtils:
    
    def __init__(self,d):
        self.d = d
        # Create a regular expression    from the dictionary keys
        #self.regex = re.compile("(%s)" % "|".join(map(re.escape, d.keys())))
        
    def multiple_replace(self, char* text): 
            """ Replace in 'text' all occurences of any key in the given
            dictionary by its corresponding value.    Returns the new string."""
            cdef int text_len = strlen(text)
            cdef char* new_text = <char *>malloc(text_len)
            strcpy(new_text, text)
            
            cdef int i 
            for i in range(len(text)):
                if self.d.has_key(chr(text[i])):
                    new_text[i] = ord(self.d[chr(text[i])])
            return new_text
            
            # For each match, look-up corresponding value in dictionary
            #return self.regex.sub(lambda mo: self.d[mo.string[mo.start():mo.end()]], text)  