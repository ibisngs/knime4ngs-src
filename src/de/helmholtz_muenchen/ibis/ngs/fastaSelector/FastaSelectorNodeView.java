package de.helmholtz_muenchen.ibis.ngs.fastaSelector;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "FastaSelector" Node.
 * This Node can be used to select multiple fasta files.
 *
 * @author Michael Kluge
 * TODO: implement the view
 */
public class FastaSelectorNodeView extends NodeView<FastaSelectorNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link FastaSelectorNodeModel})
     */
    protected FastaSelectorNodeView(final FastaSelectorNodeModel nodeModel) {
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

