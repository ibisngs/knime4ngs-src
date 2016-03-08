package de.helmholtz_muenchen.ibis.ngs.fastqc;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;


/**
 * <code>NodeDialog</code> for the "FastQC" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Max
 */
public class FastQCNodeDialog extends HTExecutorNodeDialog {

    /**
     * New pane for configuring the FastQC node.
     */
    protected FastQCNodeDialog() {

    }
}

