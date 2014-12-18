package de.helmholtz_muenchen.ibis.ngs.plotgatkgenotypeconcordance;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "PlotGATKGenotypeConcordance" Node.
 * 
 *
 * @author 
 */
public class PlotGATKGenotypeConcordanceNodeFactory 
        extends NodeFactory<PlotGATKGenotypeConcordanceNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public PlotGATKGenotypeConcordanceNodeModel createNodeModel() {
        return new PlotGATKGenotypeConcordanceNodeModel();
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
    public NodeView<PlotGATKGenotypeConcordanceNodeModel> createNodeView(final int viewIndex,
            final PlotGATKGenotypeConcordanceNodeModel nodeModel) {
        return new PlotGATKGenotypeConcordanceNodeView(nodeModel);
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
        return new PlotGATKGenotypeConcordanceNodeDialog();
    }

}

