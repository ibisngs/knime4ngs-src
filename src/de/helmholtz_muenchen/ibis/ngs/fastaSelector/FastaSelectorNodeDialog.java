package de.helmholtz_muenchen.ibis.ngs.fastaSelector;


import java.io.File;
import java.io.FilenameFilter;

import de.helmholtz_muenchen.ibis.utils.FastaFileNameFilter;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.FileSelector.FileSelectorNodeDialog;
import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;

/**
 * <code>NodeDialog</code> for the "FastaSelector" Node.
 * This Node can be used to select multiple fasta files.
 * 
 * @author Michael Kluge
 */
public class FastaSelectorNodeDialog extends FileSelectorNodeDialog {

	
	protected FastaSelectorNodeDialog() {
        super();
    }
        
	@Override
	protected String getFiletypeName() {
		return FastaSelectorNodeModel.FILTERTYPE_NAME;
	}

	@Override
	public FilenameFilter getFilenameFilter() {
		return new FastaFileNameFilter();
	}

	@Override
	public boolean isFileValid(File file) {
		// check for file ending
    	if(!new FastaFileNameFilter().accept(file.getParentFile(), file.getName()))
    		return false;
    	
    	// check, if fasta file is valid
    	if(!FileValidator.checkFastaFormat(file.getAbsolutePath()))
    		return false;
    	
    	// all checks were ok.
    	return true;
	}
}
