package de.helmholtz_muenchen.ibis.ngs.DESeq;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeView;

/**
 * <code>NodeFactory</code> for the "DESeq" Node.
 * 
 *
 * @author Michael Kluge
 */
public class DESeqNodeFactory 
        extends NodeFactory<DESeqNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public DESeqNodeModel createNodeModel() {
        return new DESeqNodeModel();
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
    public NodeView<DESeqNodeModel> createNodeView(final int viewIndex,
            final DESeqNodeModel nodeModel) {
        return new RNodeView<DESeqNodeModel>(nodeModel);
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
        return new DESeqNodeDialog();
    }

}

