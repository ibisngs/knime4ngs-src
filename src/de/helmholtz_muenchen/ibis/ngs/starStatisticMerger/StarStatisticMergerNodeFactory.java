package de.helmholtz_muenchen.ibis.ngs.starStatisticMerger;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "StarStatisticMergerNode" Node.
 * 
 *
 * @author Michael Kluge
 */
public class StarStatisticMergerNodeFactory 
        extends NodeFactory<StarStatisticMergerNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public StarStatisticMergerNodeModel createNodeModel() {
        return new StarStatisticMergerNodeModel();
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
    public NodeView<StarStatisticMergerNodeModel> createNodeView(final int viewIndex, final StarStatisticMergerNodeModel nodeModel) {
        return new StarStatisticMergerNodeView(nodeModel);
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
        return new StarStatisticMergerNodeDialog();
    }

}

