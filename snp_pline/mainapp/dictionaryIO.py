class DictionaryIO:
    def nested_dict_pairs_iterator(self, dict_obj):
        ''' This function accepts a nested dictionary as argument
            and iterate over all values of nested dictionaries
        '''
        # Iterate over all key-value pairs of dict argument
        for key, value in dict_obj.items():
            # Check if value is of dict type
            if isinstance(value, dict):
                # If value is dict then iterate over all its values
                for pair in self.nested_dict_pairs_iterator(value):
                    yield (key, *pair)
            else:
                # If value is not dict type then yield the value
                yield (key, value)
