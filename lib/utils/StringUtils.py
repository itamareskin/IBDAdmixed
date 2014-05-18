'''
Created on Dec 11, 2011

@author: itamares
'''

import re

def multiple_replace(d, text): 
        """ Replace in 'text' all occurences of any key in the given
        dictionary by its corresponding value.    Returns the new tring.""" 
    
        # Create a regular expression    from the dictionary keys
        regex = re.compile("(%s)" % "|".join(map(re.escape, d.keys())))
    
        # For each match, look-up corresponding value in dictionary
        return regex.sub(lambda mo: d[mo.string[mo.start():mo.end()]], text) 