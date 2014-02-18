package de.helmholtz_muenchen.ibis.ngs.star;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.ExecutorNode.ExecutorNodeView;

/**
 * <code>NodeView</code> for the "Star" Node.
 * STAR aligns RNA-seq reads to a reference genome using uncompressed suffix arrays. 
 * For details, please see paper: 
 * Dobin et al, Bioinformatics 2012; doi: 10.1093/bioinformatics/bts635
 *
 * @author Michael Kluge
 */
public class StarNodeView extends ExecutorNodeView<StarNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link StarNodeModel})
     */
    protected StarNodeView(final StarNodeModel nodeModel) {
        super(nodeModel);
    }
}

