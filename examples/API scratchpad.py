
q = Q3C()
q.database = "wise_index.sqlite3" # interp as SQLite3 db
q.ra_key = "ra"
q.dec_key = "dec"

q = Q3C()
q.database = <np structured table>
q.ra_key = "ra"
q.dec_key = "dec"

q.data_source = 'wise_index.sqlite3"
# + 
q.ra_key = "ra"
q.dec_key = "dec"
#or
q.data_source = <np array>
q.ra_key = 0 # index
q.dec_key = 1

q.data_source = <dict>
q.ra_key = "ra"
q.dec_key = "dec"

# or let Q3C choose
q.add_points = # list of (ra,dec) tuples
q.add_point = [ra,dec] # or (ra,dec)
# save db as pyq3c_index_tmp.sqlite3
# raise exception with provided data sources
