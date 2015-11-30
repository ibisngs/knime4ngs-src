package de.helmholtz_muenchen.ibis.ngs.vcfmerger;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "VCFMerger" Node.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class VCFMergerNodeFactory 
        extends NodeFactory<VCFMergerNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public VCFMergerNodeModel createNodeModel() {
        return new VCFMergerNodeModel(1, 1);
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
    public NodeView<VCFMergerNodeModel> createNodeView(final int viewIndex,
            final VCFMergerNodeModel nodeModel) {
        return new VCFMergerNodeView(nodeModel);
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
        return new VCFMergerNodeDialog();
    }

}

