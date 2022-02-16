suppressPackageStartupMessages({
        library(shape)
})

axis_fontsize = 0.82
small_plot_axis_fontsize = 0.715
num_lab_fontsize = 0.61
genlwd = 1

## Utility function to convert names from, say, Europe to Eur.
shorten_names = function (df) {
        names(df) = sapply(names(df),function(x) {
                gsub('Europe',     'Eur.',    x) ->.
                gsub('Africa',     'Afr.',    .) ->.
                gsub('N.Amer',     'N.Am.',   .) ->.
                gsub('S.Amer',     'S.Am.',   .) ->.
                gsub('Oceania',    'Oc.',     .) ->.
                gsub(' → ',         '→',      .) ->.
                return(.)
        })
        df
}

zoominbox = function (stub, boxmerged) {
        boxmerged = shorten_names(boxmerged)
        par(lwd = genlwd)
        if (stub == "lambsame-withmu") {
                ## Draw base box plot with grids.
                oldmgp = par('mgp')
                oldtck = par('tck')
                oldmar = par('mar')
                par(mar=oldmar*c(0.45,1.16,0.6,0.05), tck = .01, mgp=oldmgp*c(1,0.3,1))
                boxplot(boxmerged, horizontal=T, las=1, pch=20, outcex=0.6, boxwex=0.4, border = NA,
                        axes=F)
                grid(26,26)
                #title("Posterior Marginals: Model IV")
                
                ## Draw the actual main box plots, with a rectangle enclosing the transition rates.
                ## rect(xleft=-0.07, ybottom=5.5, xright=0.94, ytop=36.5, lty=2)
                segments(x0=-0.07, y0=5.5, x1=0.85,  y1=5.5, lty=3)    ## bottom segment
                segments(x0=-0.07, y0=36.5, x1=0.85, y1=36.5, lty=3)   ## top segment
                segments(x0=-0.07, y0=5.5, x1=-0.07, y1=36.5, lty=3)   ## left segment
                segments(x0=0.85, y0=5.5, x1=0.85,  y1=8., lty=3)    ## right segment lower part
                segments(x0=0.85, y0=36.5, x1=0.85,  y1=32.2, lty=3)     ## right segment upper part
                
                boxplot(boxmerged, horizontal=T, las=1, pch=20, outcex=0.6, boxwex=0.4, outline=F,
                        col="white", medlty=0, add=T,
                        xaxt   = 'n',
                        yaxt   = 'n')
                axis(2, las=1, at=seq_len(length(boxmerged)), labels=names(boxmerged),
                     cex.axis=axis_fontsize * 0.82, lwd = genlwd)
                axis(1, cex.axis=axis_fontsize * 0.82, lwd = genlwd)

                                                
                ## Draw number labels
                df = as.data.frame(boxmerged)
                colmns = colMeans(df)
                points(x=colmns, y=1:ncol(df), pch='.', cex=2.)
                text(x=colmns, y=(1:ncol(df))+0.592,
                     labels={
                         s = sprintf("%.03f",colmns)
                         s[6:length(s)] = ''
                         s
                     }, cex= num_lab_fontsize*0.9)
                                 
                ## Draw arrow
                Arrows(x0=0.25, y0=15.5, x1 = 0.52, y1 = 17.0, lwd=1.1, arr.type="triangle")
                
                ## Draw the zoomed-in plot
                par(new=T, fig=c(0.475,0.996,0.22, 0.845), tck= .007, mgp=oldmgp*c(1, 0.2, 1), mar = c(1,0.6,0.5,0.3))
                print(names(boxmerged))
                boxplot(boxmerged[6:35], horizontal=T, las=1, pch=20, outcex=0.6, boxwex=0.4, outline=F,
                                                border = "#000000FF",
                                                col    = "#FFFFFFFF",
                                                xaxt   = 'n',
                                                yaxt   = 'n',
                                                medlty = 0,
                                                axes   = F,
                                                ylim   = c(0.0, 0.85))
                axis(2, at=1:30, labels=names(boxmerged)[6:35], las=1, cex.axis=small_plot_axis_fontsize*0.7,
                     lwd=genlwd, font=2)
                axis(1, at=pretty(c(0,0.85)), labels=format(pretty(c(0,0.85),digits=2)),
                     cex.axis=small_plot_axis_fontsize*0.85, padj=-1.75, lwd=genlwd,
                     font=2)
                points(x=colmns[6:35], y=1:30, pch='.', cex=2.)
                par(tck = oldtck, mgp=oldmgp, mar=oldmar)
                ;;
        } else if (stub == "lambsameqpair-withmu") {
                ;;
        } else if (stub == "qpair-withmu") {
                ## Draw base box plot with grids.
                oldmgp = par('mgp')
                oldtck = par('tck')
                oldmar = par('mar')
                par(mar=oldmar*c(0.45,1.3,0.6,0.05), tck = .01, mgp=oldmgp*c(1,0.3,1))
                boxplot(boxmerged, horizontal=T, las=1, pch=20, outcex=0.6, boxwex=0.4, border = NA,
                        axes=F)
                grid(26,26)
                #title("Posterior Marginals: Model II")
                
                ## Draw the actual main boxes.
                rect(xleft=-0.07, ybottom=10.5, xright=0.37, ytop=26, lty=3)
                boxplot(boxmerged, horizontal=T, las=1, pch=20, outcex=0.6, boxwex=0.4, outline=F,
                                                col="white", medlty=0, add=T, cex.axis=axis_fontsize)
                                                
                ## Draw number labels and dot inside boxes
                df = as.data.frame(boxmerged)
                colmns = colMeans(df)
                points(x=colmns, y=1:ncol(df), pch='.', cex=2.)
                text(x=colmns, y=(1:ncol(df))+0.575,
                    labels={
                               s = sprintf("%.03f",colmns)
                               s[11:length(s)] = ''
                               s
                           }, cex=num_lab_fontsize)
                                 
                ## Draw arrow
                Arrows(x0=2/8, y0=11.6, x1 = 0.42, y1 = 13.0, lwd=1.1, arr.type="triangle")
                
                ## Draw the zoomed-in plot
                par(new=T, fig=c(0.5,1.0,0.43, 0.925), tck= .007, mgp=oldmgp*c(0.1, 0.1, 0.1), mar = c(1,0.6,0.5,0.3))
                boxplot(boxmerged[11:25], horizontal=T, las=1, pch=20, outcex=0.6, boxwex=0.4, outline=F,
                        border = "#000000FF",
                        col    = "#FFFFFFFF",
                        xaxt   = 'n',
                        yaxt   = 'n',
                        medlty = 0,
                        axes   = F,
                        ylim   = c(0.0, 0.35))
                axis(2, at=1:15, labels=names(boxmerged)[11:25], las=1, cex.axis=small_plot_axis_fontsize)
                axis(1, at=pretty(c(0,0.35)), labels=format(pretty(c(0,0.35),digits=2)),
                     cex.axis=small_plot_axis_fontsize, padj=-.8)
                points(x=colmns[11:25], y=1:15, pch='.', cex=2.)
                par(tck = oldtck, mgp=oldmgp, mar=oldmar)
        } else {
                ;;
        }
}
