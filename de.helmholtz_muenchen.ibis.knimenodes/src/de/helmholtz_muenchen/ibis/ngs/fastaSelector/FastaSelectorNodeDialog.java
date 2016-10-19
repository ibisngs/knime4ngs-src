/**
 *  Copyright (C) 2016 the Knime4NGS contributors.
 *  Website: http://ibisngs.github.io/knime4ngs
 *  
 *  This file is part of the KNIME4NGS KNIME extension.
 *  
 *  The KNIME4NGS extension is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package de.helmholtz_muenchen.ibis.ngs.fastaSelector;


import java.io.File;

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
	public boolean isFileValid(File file) {
		// check for file ending
    	if(!this.getFilenameFilter().accept(file.getParentFile(), file.getName()))
    		return false;
    	
    	// check, if fasta file is valid
    	if(!FileValidator.checkFastaFormat(file.getAbsolutePath()))
    		return false;
    	
    	// all checks were ok.
    	return true;
	}

	@Override
	public String getFilenameEndRegex() {
		return FastaSelectorNodeModel.FILENAME_END_REGEX;
	}
}
