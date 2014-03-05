package de.helmholtz_muenchen.ibis.ngs.samSelector;


import java.io.File;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.FileSelector.FileSelectorNodeDialog;
import de.helmholtz_muenchen.ibis.utils.ngs.FileValidator;

/**
 * <code>NodeDialog</code> for the "FastaSelector" Node.
 * This Node can be used to select multiple fasta files.
 * 
 * @author Michael Kluge
 */
public class SamSelectorNodeDialog extends FileSelectorNodeDialog {

	
	protected SamSelectorNodeDialog() {
        super();
    }
        
	@Override
	protected String getFiletypeName() {
		return SamSelectorNodeModel.FILTERTYPE_NAME;
	}

	@Override
	public boolean isFileValid(File file) {
		// check for file ending
    	if(!this.getFilenameFilter().accept(file.getParentFile(), file.getName()))
    		return false;
    	
    	// TODO check, if fasta file is valid
    	//if(!FileValidator. .checkFastaFormat(file.getAbsolutePath()))
    	//	return false;
    	
    	// all checks were ok.
    	return true;
	}

	@Override
	public String getFilenameEndRegex() {
		return SamSelectorNodeModel.FILENAME_END_REGEX;
	}
}