package de.helmholtz_muenchen.ibis.ngs.DESeq;

import org.knime.core.node.NodeView;

/**
 * <code>NodeView</code> for the "DESeq" Node.
 * 
 *
 * @author Michael Kluge
 */
public class DESeqNodeView extends NodeView<DESeqNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link DESeqNodeModel})
     */
    protected DESeqNodeView(final DESeqNodeModel nodeModel) {
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

