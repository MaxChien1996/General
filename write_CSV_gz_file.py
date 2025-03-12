import csv
import gzip



def write_CSV_gz_file(dict, name_of_1st_column, output_CSV_gz):
    """Write a dictionary to a CSV.gz file"""
    """By dynamically detecting the fieldnames (keys of the dictionary)"""
    """However, due to the nature of Python dictionary, it needs the first fieldname by manually input."""
    
    if not dict:
        
        raise ValueError("The input dictionary is empty. Cannot write to CSV.")
    
    
    
    # Dynamically extract field names from the dictionary keys, but need to manually input the first fieldname
    first_key = next(iter(dict))
        # iter() returns an iterator over the keys of the dictionary
        # next() retrieves the first key from the iterator
    fieldnames = [name_of_1st_column] + list(dict[first_key].keys())
        # due to the lack of first column name in the dictionary, we add it manually
    
    
    
    with gzip.open(output_CSV_gz, 'wt', newline='') as CSV_file:
        
        writer = csv.DictWriter(CSV_file, fieldnames = fieldnames)
        writer.writeheader()
        
        for key, item in dict.items():
            
            row = {name_of_1st_column: key}
            row.update(item)  # This adds all key-value pairs from `item` into `row`
            writer.writerow(row)