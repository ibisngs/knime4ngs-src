package de.helmholtz_muenchen.ibis.ngs.kggseq;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "KGGSeq" Node.
 * 
 *
 * @author 
 */
public class KGGSeqNodeFactory 
        extends NodeFactory<KGGSeqNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public KGGSeqNodeModel createNodeModel() {
        return new KGGSeqNodeModel();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNrNodeViews() {
        return 1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeView<KGGSeqNodeModel> createNodeView(final int viewIndex,
            final KGGSeqNodeModel nodeModel) {
        return new KGGSeqNodeView(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean hasDialog() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeDialogPane createNodeDialogPane() {
        return new KGGSeqNodeDialog();
    }

}

