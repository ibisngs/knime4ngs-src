package de.helmholtz_muenchen.ibis.ngs.samSelector;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "FastaSelector" Node.
 * This Node can be used to select multiple fasta files.
 *
 * @author Michael Kluge
 * TODO: implement the view that it displays the selected fasta files which where written to the output
 */
public class SamSelectorNodeView extends NodeView<SamSelectorNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link SamSelectorNodeModel})
     */
    protected SamSelectorNodeView(final SamSelectorNodeModel nodeModel) {
        super(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void modelChanged() { 
    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onClose() {
    	
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected void onOpen() {

    }
}

