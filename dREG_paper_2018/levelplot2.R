
## version using contourLines, and hopefully works for missing matrix
## entries as well

panel.levelplot <-
    function(x, y, z,
             subscripts,
             at = pretty(z),
             shrink,
             labels = FALSE,
             label.style = c("mixed", "flat", "align"),
             contour = FALSE,
             region = TRUE,
             col = add.line$col,
             lty = add.line$lty,
             lwd = add.line$lwd,
             border = "transparent",
             border.lty = 1,
             border.lwd = 0.1,
             ...,
             col.regions = regions$col,
             alpha.regions = regions$alpha,
             identifier = "levelplot")
{
    if (length(subscripts) == 0) return()
    regions <- trellis.par.get("regions")
    label.style <- match.arg(label.style)
    x.is.factor <- is.factor(x)
    y.is.factor <- is.factor(y)
    x <- as.numeric(x)
    y <- as.numeric(y)
    z <- as.numeric(z)

##     numcol <- length(at) - 1
##     numcol.r <- length(col.regions)
##     col.regions <-
##         if (numcol.r <= numcol)
##             rep(col.regions, length.out = numcol)
##         else
##             col.regions[round(seq(1, numcol.r, length.out = numcol))]
##     zcol <- cut(z, at, include.lowest = TRUE, labels = FALSE)

    zcol <-
        if (region) level.colors(z, at, col.regions, colors = TRUE) else "transparent"

    x <- x[subscripts]
    y <- y[subscripts]
    minXwid <- if (length(unique(x)) > 1) min(diff(sort(unique(x)))) else 1
    minYwid <- if (length(unique(y)) > 1) min(diff(sort(unique(y)))) else 1
    fullZrange <- range(as.numeric(z), finite = TRUE) # for shrinking
    z <- z[subscripts]
    if (region) zcol <- zcol[subscripts]

    if (lattice:::hasGroupNumber()) group <- list(...)$group.number  else       group <- 0

    ## Do we need a zlim-like argument ?

    shrinkx <- c(1, 1)
    shrinky <- c(1, 1)
    if (!missing(shrink)) {
        if (is.numeric(shrink)) {
            shrinkx <- rep(shrink, length.out = 2)
            shrinky <- rep(shrink, length.out = 2)
        }
        else if (is.list(shrink)) {
            shrinkx <- rep(shrink[[1]], length.out = 2)
            shrinky <- rep(shrink[[1]], length.out = 2)
            if ("x" %in% names(shrink)) shrinkx <- rep(shrink$x, length.out = 2)
            if ("y" %in% names(shrink)) shrinky <- rep(shrink$y, length.out = 2)
        }
        else warning("Invalid 'shrink' parameter ignored")
    }

    scaleWidth <- function(z, min = .8, max = .8, zl = range(z, finite = TRUE)) {
        if (diff(zl) == 0) rep(.5 * (min + max), length(z))
        else min + (max - min) * (z - zl[1]) / diff(zl)
    }

    if (x.is.factor)
    {
        ## unique values
        ux <- sort(unique(x[!is.na(x)]))
        ## dimension of rectangles
        lx <- rep(1, length(ux))
        ## centers of rectangles
        cx <- ux
    }
    else
    {
        ## sorted unique values of x
        ux <- sort(unique(x[!is.na(x)]))
        ## actual box boundaries (x axis)
        bx <-
            if (length(ux) > 1)
                c(3 * ux[1] - ux[2],
                  ux[-length(ux)] + ux[-1],
                  3 * ux[length(ux)] - ux[length(ux)-1]) / 2
            else
                ux + c(-.5, .5) * minXwid
        ## dimension of rectangles
        lx <- diff(bx)
        ## centers of rectangles
        cx <- (bx[-1] + bx[-length(bx)])/2
    }

	## same things for y
    if (y.is.factor)
    {
        ## unique values
        uy <- sort(unique(y[!is.na(y)]))
        ## dimension of rectangles
        ly <- rep(1, length(uy))
        ## centers of rectangles
        cy <- uy
    }
    else
    {
        uy <- sort(unique(y[!is.na(y)]))
        by <-
            if (length(uy) > 1)
                c(3 * uy[1] - uy[2],
                  uy[-length(uy)] + uy[-1],
                  3 * uy[length(uy)] - uy[length(uy)-1]) / 2
            else
                uy + c(-.5, .5) * minYwid
        ly <- diff(by)
        cy <- (by[-1] + by[-length(by)])/2
    }

    idx <- match(x, ux)
    idy <- match(y, uy)

    if (region)
	{
		for(i in 1:NROW(idx))
		{
			if( idy[i]>=idx[i] )
		        grid.rect(x = cx[idx[i]],
	                  y = cy[idy[i]],
	                  width = lx[idx[i]] * scaleWidth(z[i], shrinkx[1], shrinkx[2], fullZrange),
	                  height = ly[idy[i]] * scaleWidth(z[i], shrinky[1], shrinky[2], fullZrange),
	                  default.units = "native",
	                  name = trellis.grobname(paste(identifier, "rect", sep="."), type = "panel", group = group),
	                  gp =
	                  gpar(fill = zcol[i],
	                       col = border,
	                       lwd = border.lwd,
	                       lty = border.lty,
	                       alpha = alpha.regions))
	        else
	        {
                width = lx[idx[i]] * scaleWidth(z[i], shrinkx[1], shrinkx[2], fullZrange);
	            height = ly[idy[i]] * scaleWidth(z[i], shrinky[1], shrinky[2], fullZrange);

		        grid.rect(x = cx[idx[i]],
	                  y = cy[idy[i]],
	                  width = width,
	                  height = height,
	                  default.units = "native",
	                  name = trellis.grobname(paste(identifier, "rect", sep="."), type = "panel", group = group),
	                  gp =
	                  gpar(fill = NA,
	                       col = zcol[1],
	                       lwd = 1,
	                       alpha = alpha.regions));

				dr <- range(g_pairs[!is.infinite(g_pairs)]);
				dx <- (g_pairs[,idx[i]]-dr[1])/(dr[2]-dr[1]);
				dy <- (g_pairs[,idy[i]]-dr[1])/(dr[2]-dr[1]);

				if(!exists("g_pairs_col")) g_pairs_col <<- zcol[1];
				if(!exists("g_pairs_size")) g_pairs_col <<- 0.01;

	        	grid.points(x = cx[idx[i]] - 0.4*width + dx*(width*0.8),
				            y = cy[idy[i]] - 0.4*height + dy*(height*0.8),
				            pch = 21,
				            size = unit(g_pairs_size, "char"),
				            gp = gpar( col = g_pairs_col ) );
			}
	 	}
	}

    if (contour)
    {
        ## calculate aspect ratio of panel to use in calculating label alignment
        cpl <- current.panel.limits(unit="cm")
        asp <- diff(cpl$ylim) / diff(cpl$xlim)

        ## Processing the labels argument
        if (is.logical(labels) && !labels) labels <- NULL
        else
        {
            if (is.characterOrExpression(labels)) labels <- list(labels = labels)
            text <- trellis.par.get("add.text")
            tmp <- list(col = text$col,
                        alpha = text$alpha,
                        cex = text$cex,
                        fontfamily = text$fontfamily,
                        fontface = text$fontface,
                        font = text$font)
            labels <-
                if (is.list(labels)) updateList(tmp, labels)
                else tmp
            if (!is.characterOrExpression(labels$labels)) # NULL/TRUE
                labels$labels <- format(at, trim = TRUE)
        }

        add.line <- trellis.par.get("add.line")

        ## convert z into a matrix, with NA entries for those
        ## 'missing' from data frame. There's scope for ambiguity
        ## here, which can be avoided by the user.

        m <- matrix(NA_real_, nrow = length(ux), ncol = length(uy))
        m[(idy - 1) * length(ux) + idx ] <- z

        clines <-
            contourLines(x = ux, y = uy, z = m,
                         nlevels = length(at), ## necessary ?
                         levels = at)

        ccount <- 0

        for (val in clines) {

            ccount <- ccount + 1

            ## each val looks like:

            ## $ :List of 3
            ##  ..$ level: num 170
            ##  ..$ x    : num [1:21] 0.535 0.534 0.534 0.534 0.535 ...
            ##  ..$ y    : num [1:21] 0.398 0.400 0.417 0.433 0.434 ...

            ## we don't know how to leave gap in lines for labels.

            llines(val, ## hopefully $levels won't matter
                   col = col, lty = lty, lwd = lwd,
                   identifier = paste(identifier, "line", ccount,
                     sep = "."))

            ## if too small, don't add label. How small is small ?
            ## Should depend on resolution. How ?

            if (length(val$x) > 5)
            {
                if (!is.null(labels))
                {
                    slopes <- diff(val$y) / diff(val$x)
                    ## slopes[is.na(slopes)] <- 0

                    if (label.style == "flat")
                    {
                        ## draw label at 'flattest' position along contour

                        textloc <- which.min(abs(slopes))
                        rotangle <- 0
                    }
                    else if (label.style == "align")
                    {

                        ## draw label at 'deepest' position along
                        ## contour, depth being min distance to either
                        ## of the four edges, scaled appropriately

                        rx <- range(ux)
                        ry <- range(uy)
                        depth <- pmin(pmin(val$x - rx[1], rx[2] - val$x) / diff(rx),
                                      pmin(val$y - ry[1], ry[2] - val$y) / diff(ry))
                        textloc <- min(which.max(depth), length(slopes))
                                        # slopes has one less entry,
                                        # and textloc indexes slopes

                        rotangle <- atan(asp * slopes[textloc] * diff(rx) / diff(ry)) * 180 / base::pi
                    }
                    else if (label.style == "mixed")
                    {

                        ## mix both. align for contours whose flattest
                        ## portion is too close to edge

                        rx <- range(ux)
                        ry <- range(uy)
                        depth <- pmin(pmin(val$x - rx[1], rx[2] - val$x) / diff(rx),
                                      pmin(val$y - ry[1], ry[2] - val$y) / diff(ry))
                        textloc <- which.min(abs(slopes))
                        rotangle <- 0

                        if (depth[textloc] < .05 ) {
                            textloc <- min(which.max(depth), length(slopes))
                            rotangle <- atan(asp * slopes[textloc] * diff(rx) / diff(ry)) * 180 / base::pi
                        }
                    }
                    else stop("Invalid label.style")

                    i <- match(val$level, at)

                    ltext(labels$labels[i],
                          adj = c(.5, 0),
                          srt = rotangle,
                          col = labels$col,
                          alpha = labels$alpha,
                          cex = labels$cex,
                          font = labels$font,
                          fontfamily = labels$fontfamily,
                          fontface = labels$fontface,
                          x = .5 * (val$x[textloc]+val$x[textloc + 1]),
                          y = .5 * (val$y[textloc]+val$y[textloc + 1]),
                          identifier = paste(identifier, "label", ccount,
                            sep = "."))
                }
            }
        }
    }
}
