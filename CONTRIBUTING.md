## Code Documentation

Our SONATA project uses the Syntax and best practices for docstrings to be used with the numpydoc extension for [Sphinx](http://sphinx-doc.org/). Allowing our tool to produce well-formatted reference guides. Our docstring standard uses [re-structured text (reST)](http://docutils.sourceforge.net/rst.html) syntax and is rendered using [Sphinx](http://sphinx.pocoo.org/) (a pre-processor that understands the particular documentation style we are using). Keep in mind that the length of docstring lines should be kept to 75 characters to facilitate reading the docstrings in text terminals.

The complete Styleguide and can be found here: https://numpydoc.readthedocs.io/en/latest/format.html with an of snippets with [example.py](example.py)

The following example of the CBM method cbm_run_vabs illustrates the idea:

```python
        def cbm_run_vabs(self, jobid=None, rm_vabfiles=True, ramdisk=False):
        '''CBM method to run the solver VABS (Variational Asymptotic Beam 
        Sectional Analysis). Note that this method is designed to work if 
        VABSIII is set in the PATH variable. For Users at the TUM-HT please load 
        the vabs module beforehand.
                
        Parameters
        ----------
        jobid : string, optional
                assign a unique ID for the job. If no jobid is assigned the 
                isoformat of datetime with microseconds is used
        rm_vabfiles : bool, optional
                removes VABS files after the calculation is completed and 
                the results are stored.
        ramdisk : bool, optional, 
            Instead of storing the writing and reading the vabs job directory, 
            the ramdisk "/tmpfs/username" is used. This options is currently 
            designed for linux users make sure to mount it beforehand with to 
            assign 200MB of Memory to the virtual drive.
            >>> sudo mount -t tmpfs -o size=200M none /tmpfs/username
            
        Returns
        ----------
        None : everything is stored within the CBM instance
        
        Examples
        ----------
        >>> job.cbm_run_vabs(rm_vabfiles=True, ramdisk=True)

        '''
```