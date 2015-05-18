package de.helmholtz_muenchen.ibis.ngs.vcfnormalizer;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "VCFNormalizer" Node.
 * 
 *
 * @author tim.jeske
 */
public class VCFNormalizerNodeFactory 
        extends NodeFactory<VCFNormalizerNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public VCFNormalizerNodeModel createNodeModel() {
        return new VCFNormalizerNodeModel();
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
    public NodeView<VCFNormalizerNodeModel> createNodeView(final int viewIndex,
            final VCFNormalizerNodeModel nodeModel) {
        return new VCFNormalizerNodeView(nodeModel);
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
        return new VCFNormalizerNodeDialog();
    }

}

