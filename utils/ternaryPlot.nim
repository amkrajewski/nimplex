import nimplex
import plotting
import arraymancer/Tensor
import pixie
import std/strutils
import std/math
import std/times

let t0 = cpuTime()

# User-set Constants
const
    nDiv: int = 24
    filename: string = "output.png"
    height: int = 4000
    scaleChroma: float = 2.0
    fontName: string = "MM" # "NASA", "IBM", "DM", "MM"
    indexOverlay: bool = false
    compIsMain: bool = false
    marker1: string = "Undesirable" # Yellow
    markerOverlay1: bool = true
    marker2: string = "Infeasible"  # Red
    markerOverlay2: bool = true
    propertyOverlay: bool = false
    pathPointsOverlay: bool = true
    pathType: string = "highvis" # "line" or "boldline" or "highvis"
    propetyColoringStyle: string = "hotcold"
    elementList: seq[string] = @["Fe", "Cr", "Ni", "Ti", "Cu"] # e.g., @["Fe", "Cr", "Ni", "Ti", "Cu"] or @["Ti", "Mo", "Cu"]
    labels: seq[string] = @["Ti64", "SS316L", "Monel400"]
    elementalCompositions: seq[seq[float]] = @[@[0.0, 0.0, 0.0, 100.0, 0.0],  @[68.0, 18.0, 14.0, 0.0, 0.0],  @[2.0, 0.0, 66.0, 0.0, 32.0]]
    # elementalCompositions: seq[seq[float]] = @[@[1,0,0], @[0,1,0], @[0,0,1]]
    # @[@[20.0, 40.0, 0.0, 0.0, 0.0], @[0.5, 0.0, 0.0, 0.0, 0.5], @[0.1, 0.0, 0.5, 0.5, 0.1]])
    axisType: string = "vulgarFractions" # "decimal" or "vulgarFractions" or "quanta"
    pathPoints: seq[int] = @[0, 94, 95, 173, 174, 222, 223, 283, 286, 267, 268, 245, 246, 24]
    scaling: float = height/2000

# Nimplex Processing
let 
    grid: Tensor[float] = nimplex.simplex_grid_fractional(3, nDiv)
    elementalEmbedding: Tensor[float] = nimplex.attainable2elemental(grid, elementalCompositions) 
    pureCompIDX = nimplex.pure_component_indexes(3, nDiv)
var gridPoints = plotting.simplex2cartesian(grid) * (1000 * scaling)
gridPoints[_, 0] = gridPoints[_, 0] +. (1250 * scaling)
gridPoints[_, 1] = -gridPoints[_, 1] +. (1250 * scaling)

let
    gpl = gridPoints.toSeq2D
    (pp1, pp2, pp3) = (gpl[pureCompIDX[0]], gpl[pureCompIDX[1]], gpl[pureCompIDX[2]])

# Derived Constants
const
    elementalEmbeddingLen: int = elementalCompositions[0].len
    longLabel: bool = (len(labels[1]) > 2 or len(labels[2]) > 2)
    width: int = ((height.float)*(5/4)).int
    sideWidth: float = (60*scaling)
    luminance: float = if compIsMain: 0.55 else: 0.45
    alpha: uint8 = if compIsMain: 255 else: 128
    # Color schemes
    colorSchemeOKlab2c: seq[seq[float]] = @[
        @[luminance,  0.250*scaleChroma,  0.000],
        @[luminance, -0.250*scaleChroma,  0.000]]
    colorSchemeOKlab3c: seq[seq[float]] = @[
        @[luminance,  0.217*scaleChroma,  0.125*scaleChroma],
        @[luminance, -0.065*scaleChroma, -0.241*scaleChroma],
        @[luminance, -0.177*scaleChroma,  0.177*scaleChroma]]
    colorSchemeOKlab4c: seq[seq[float]] = @[
        @[luminance,  0.250*scaleChroma,  0.000],
        @[luminance, -0.250*scaleChroma,  0.000],
        @[luminance,  0.000, -0.250*scaleChroma],
        @[luminance,  0.000,  0.250*scaleChroma]]
    colorSchemeOKlab5c: seq[seq[float]] = @[
        @[luminance,  0.217*scaleChroma,  0.125*scaleChroma],
        @[luminance,  0.177*scaleChroma, -0.177*scaleChroma],
        @[luminance, -0.125*scaleChroma, -0.217*scaleChroma],
        @[luminance, -0.241*scaleChroma,  0.065*scaleChroma],
        @[luminance,  0.000*scaleChroma,  0.250*scaleChroma]]
    propertyColoringOKlabHotCold: seq[seq[float]] = @[
        @[0.700, 0.217*scaleChroma, 0.125*scaleChroma],
        @[0.7, -0.217*scaleChroma, -0.125*scaleChroma]]
    propertyColoringOKlab01: seq[seq[float]] = @[
        @[0.400, 0.217*scaleChroma, 0.125*scaleChroma],
        @[1.0, 0.0, 0.0]]

# *** Color Schemes ***
const colorSchemeOKlab: seq[seq[float]] = static:
    case elementalEmbeddingLen:
        of 2: colorSchemeOKlab2c
        of 3: colorSchemeOKlab3c
        of 4: colorSchemeOKlab4c
        of 5: colorSchemeOKlab5c
        else: raise newException(ValueError, "Unsupported number of components")

const propertyColoringOKlab: seq[seq[float]] = static:
    case propetyColoringStyle:
        of "hotcold": propertyColoringOKlabHotCold
        of "01": propertyColoringOKlab01
        else: raise newException(ValueError, "Unsupported property coloring style")

# *** Fonts ***
const 
    fontMain: string = static:
        case fontName
            of "NASA": "fonts/nasalization-rg.otf"
            of "IBM": "fonts/IBMPlexMono-Medium.ttf"
            of "DM": "fonts/DMMono-Medium.ttf"
            of "MM": "fonts/MartianMonoSemiCondensed-Medium.otf"
            else: raise newException(ValueError, "Unknown font name")
    fontSupport: string = static:
        case fontName
            of "NASA": "fonts/nasalization-rg.otf"
            of "IBM": "fonts/IBMPlexMono-Regular.ttf"
            of "DM": "fonts/DMMono-Regular.ttf"
            of "MM": "fonts/MartianMonoCondensed-Regular.otf"
            else: raise newException(ValueError, "Unknown font name")

# *** Line Widths ***
const distance: float = static:
    if nDiv < 200:
        (990/(nDiv+1)-1)*scaling
    else:
        (990/(nDiv+1)-0.7)*scaling

const (thinLine, thickLine): (float, float) = static:
    if nDiv <= 36:
        (scaling * 10 / (nDiv/12), scaling * 20 / (nDiv/12))
    elif nDiv <= 48:
        (scaling * 3, scaling * 5)
    elif nDiv < 200:
        (scaling * 1, scaling * 3) 
    else:
        (scaling * 0.5, scaling * 2) 


# *** Helper Functions ***
func clamp(x: float, minVal: float, maxVal: float): float {.inline.} =
        return max(minVal, min(maxVal, x))

func oklabToRGB(c: tuple[L, a, b: float]): tuple[r, g, b: uint8] {.inline.} =

    func linearToGamma(c: float): float {.inline.} =
        if c >= 0.0031308: 
            return 1.055 * pow(c, 1 / 2.4) - 0.055
        return 12.92 * c

    let
        l = (c.L + c.a * +0.3963377774 + c.b * +0.2158037573) ^ 3
        m = (c.L + c.a * -0.1055613458 + c.b * -0.0638541728) ^ 3
        s = (c.L + c.a * -0.0894841775 + c.b * -1.2914855480) ^ 3
        r = (255 * linearToGamma(
            l * +4.0767416621 + m * -3.3077115913 + s * +0.2309699292)
            ).clamp(0, 255).round.uint8
        g = (255 * linearToGamma(
            l * -1.2684380046 + m * +2.6097574011 + s * -0.3413193965)
            ).clamp(0, 255).round.uint8
        b = (255 * linearToGamma(
            l * -0.0041960863 + m * -0.7034186147 + s * +1.7076147010)
            ).clamp(0, 255).round.uint8

    return (r, g, b)

# Dummy data (!!!) for plot design
var 
    propertyField: Tensor[float] = newTensor[float]([grid.shape[0]])
    feasibilityField1: Tensor[bool] = newTensor[bool]([grid.shape[0]])
    feasibilityField2: Tensor[bool] = newTensor[bool]([grid.shape[0]])
    #propertyFieldNormalized: Tensor[float] = newTensor[float]([grid.shape[0]])
    elementalColoring: seq[ColorRGBA]
    propertyColoring: seq[ColorRGBA]

for i in 0..<grid.shape[0]:
    propertyField[i] = grid[i, 0] * 0.653 + grid[i, 1] * 0.587 + grid[i, 2] * 0.114 + 0.3 * (grid[i, 0] - grid[i, 1])^2 - 0.1 * (grid[i, 0] - grid[i, 1])^3

    if (grid[i, 0] + grid[i, 1]) > 0.6 and propertyField[i]<0.6:
        feasibilityField1[i] = true
    else:
        feasibilityField1[i] = false

    if (grid[i, 1] > 0.2 and propertyField[i]<0.55) or (grid[i, 0] > 0.8):
        feasibilityField2[i] = true
    else:
        feasibilityField2[i] = false

let
    minVal = propertyField.min
    maxVal = propertyField.max
    rangeVal = maxVal - minVal
    propertyFieldNormalized: Tensor[float] = (propertyField -. minVal) /. rangeVal

# Color mixing in the OKlab space for continuous coloring. Convert to RGB for display.
for i in 0..<grid.shape[0]:
    var oklSeq: seq[float] = newSeq[float](3)
    for j in 0..<elementalEmbeddingLen:
        for k in 0..<3:
            oklSeq[k] += elementalEmbedding[i, j] * colorSchemeOKlab[j][k].float
    let rgbSeq = oklabToRGB((oklSeq[0], oklSeq[1], oklSeq[2]))
    elementalColoring.add(rgba(rgbSeq[0].uint8, rgbSeq[1].uint8, rgbSeq[2].uint8, alpha))

func prop2rgb(prop: float, propertyColoringOKlab: seq[seq[float]]): ColorRGBA =
    var oklSeq: seq[float] = newSeq[float](3)
    for k in 0..<3:
        oklSeq[k] = prop * propertyColoringOKlab[0][k] + (1 - prop) * propertyColoringOKlab[1][k]
    let rgbSeq = oklabToRGB((oklSeq[0], oklSeq[1], oklSeq[2]))
    return rgba(rgbSeq[0].uint8, rgbSeq[1].uint8, rgbSeq[2].uint8, 255.uint8)

# Convert property field to color
for i in 0..<grid.shape[0]:
    propertyColoring.add(prop2rgb(propertyFieldNormalized[i], propertyColoringOKlab))
        

# **************** Drawing Functions ****************
    
proc fillWhite(image: Image) =
    image.fill(rgba(255, 255, 255, 255))

# *** BACKGROUND ***
# Thick line in the background
proc drawBackground(image: Image) =
    let ctx = newContext(image)
    ctx.fillStyle = "#140700"
    ctx.strokeStyle = "#140700"
    ctx.lineWidth = thickLine

    ctx.strokeSegment(segment(vec2(pp1[0], pp1[1]), vec2(pp2[0], pp2[1])))
    ctx.strokeSegment(segment(vec2(pp2[0], pp2[1]), vec2(pp3[0], pp3[1])))
    ctx.strokeSegment(segment(vec2(pp3[0], pp3[1]), vec2(pp1[0], pp1[1])))


# *** MIDGROUND ***
# Main hexes and points for later steps
proc drawCompositonHexes(
        image: Image, 
        elementalColoring: seq[ColorRGBA]
        ): void =
    for i in 0..<gpl.len:
        # Main hexes filling the space with gaps
        let pathHex = newPath()
        pathHex.polygon(vec2(gpl[i][0], gpl[i][1]), distance, sides = 6)
        image.fillPath(pathHex, elementalColoring[i])

# Property overlay
proc drawPropertyHexes(
        image: Image, 
        propertyColoring: seq[ColorRGBA]
        ): void =
    # Black rim
    let pathHex = newPath()
    for i in 0..<gpl.len:
        pathHex.polygon(vec2(gpl[i][0], gpl[i][1]), distance*3.1/4, sides = 6)
    image.fillPath(pathHex, rgba(0, 0, 0, 255))
    # Property color fill
    for i in 0..<gpl.len:
        let pathHex = newPath()
        pathHex.polygon(vec2(gpl[i][0], gpl[i][1]), distance*3/4, sides = 6)
        image.fillPath(pathHex, propertyColoring[i])

# Designed path overlay
proc drawDesignedPath(
        image: Image,
        pathPoints: seq[int],
        ): void =
    let ctx = newContext(image)
    if pathType == "boldline" or pathType == "highvis":
        if pathType == "highvis":
            ctx.strokeStyle = rgba(20, 0, 180, 255)
        else:
            ctx.strokeStyle = rgba(0, 180, 0, 255)
        for i in 0..<pathPoints.len-1:
            ctx.lineWidth = distance
            ctx.strokeSegment(segment(
                vec2(gpl[pathPoints[i]][0], gpl[pathPoints[i]][1]), 
                vec2(gpl[pathPoints[i+1]][0], gpl[pathPoints[i+1]][1])))
            ctx.lineWidth = distance*0.866
            ctx.strokePolygon(vec2(gpl[pathPoints[i]][0], gpl[pathPoints[i]][1]), 0.01, sides = 6)
    ctx.strokeStyle = rgba(0, 180, 0, 255)
    ctx.lineWidth = thinLine*2
    for i in 0..<pathPoints.len-1:
        ctx.strokeSegment(segment(
            vec2(gpl[pathPoints[i]][0], gpl[pathPoints[i]][1]), 
            vec2(gpl[pathPoints[i+1]][0], gpl[pathPoints[i+1]][1])))
    ctx.strokeStyle = rgba(0, 220, 0, 255)
    ctx.lineWidth = thinLine
    for i in 0..<pathPoints.len-1:
        ctx.strokeSegment(segment(
            vec2(gpl[pathPoints[i]][0], gpl[pathPoints[i]][1]), 
            vec2(gpl[pathPoints[i+1]][0], gpl[pathPoints[i+1]][1])))

# Feasibility overlay 1 (Yellow)
proc drawMarkers1(
        image: Image, 
        markerField1: Tensor[bool]
        ): void =
    let backgroundF1 = newPath()
    for i in 0..<gpl.len:
        if markerField1[i]:
            backgroundF1.circle(gpl[i][0], gpl[i][1], distance*1.05/2)
    image.fillPath(backgroundF1, rgba(0, 0, 0, 255))
    let pathF1 = newPath()
    for i in 0..<gpl.len:
        if markerField1[i]:
            pathF1.circle(gpl[i][0], gpl[i][1], distance*1/2)
    image.fillPath(pathF1, rgba(200, 200, 0, 255))

# Feasibility overlay 2 (Red)
proc drawMarkers2(
        image: Image, 
        markerField1: Tensor[bool],
        markerField2: Tensor[bool]
        ): void =
    let backgroundF2 = newPath()
    for i in 0..<gpl.len:
        if markerField2[i] and not (markerField1[i] and markerOverlay1):
            backgroundF2.circle(gpl[i][0], gpl[i][1], distance*0.75/2)
    image.fillPath(backgroundF2, rgba(0, 0, 0, 255))
    let pathF2 = newPath()
    for i in 0..<gpl.len:
        if markerField2[i]:
            pathF2.circle(gpl[i][0], gpl[i][1], distance*0.7/2)
            image.fillPath(pathF2, rgba(255, 0, 0, 255))

proc drawMarkers2(
        image: Image,
        markerField2: Tensor[bool]
        ): void =
    let markerField1 = newTensor[bool]([grid.shape[0]]) # Initialized as false by default
    drawMarkers2(image, markerField1, markerField2)

# *** FOREGROUND ***

# Thin side line in the foreground
proc drawForeground(image: Image): void =
    let ctx = newContext(image)
    ctx.strokeStyle = "#140700"
    ctx.lineWidth = thinLine
    ctx.strokeSegment(segment(vec2(pp1[0], pp1[1]), vec2(pp2[0], pp2[1])))
    ctx.strokeSegment(segment(vec2(pp2[0], pp2[1]), vec2(pp3[0], pp3[1])))
    ctx.strokeSegment(segment(vec2(pp3[0], pp3[1]), vec2(pp1[0], pp1[1])))

    let centerPoints = newPath()
    for p in gpl:
        if nDiv < 48:
            # Tiny hexes if they can be seen
            centerPoints.polygon(vec2(p[0], p[1]), thinLine, 6)
        else:
            # Otherwise, just points
            centerPoints.circle(p[0], p[1], thinLine)
    image.fillPath(centerPoints, rgba(0, 10, 10, 255))

    # White side line
    ctx.strokeStyle = rgba(255, 255, 255, 200)
    ctx.lineWidth = 2*sideWidth-thickLine/2

    ctx.strokeSegment(
        segment(
            vec2(pp1[0]-80*scaling, pp1[1]-2*sideWidth-140*scaling), 
            vec2(pp2[0]+80*scaling+2*sideWidth*0.58, pp2[1]+140*scaling)))
    ctx.strokeSegment(
        segment(
            vec2(pp2[0]+2*sideWidth, pp2[1]+2*sideWidth/2), 
            vec2(pp3[0]-2*sideWidth, pp3[1]+2*sideWidth/2)))
    ctx.strokeSegment(
        segment(
            vec2(pp3[0]-80*scaling-2*sideWidth*0.58, pp3[1]+140*scaling), 
            vec2(pp1[0]+80*scaling, pp1[1]-2*sideWidth-140*scaling)))

# ********* Axis Labels *********

proc drawAxisLabels(image: Image) =
    var font: Font = readFont(fontMain)
    font.size = sideWidth*1.5
    font.paint.color = color(0.7, 0, 0.1)
    # Top label is always top-center aligned and needs no horizontal tuning.
    image.fillText(
        font,
        labels[0], 
        translate(vec2(pp1[0], pp1[1]-sideWidth*2)),
        halign = CenterAlign
        )

    if longLabel:
        # Longer labels need to live under the axis
        image.fillText(font, labels[1], translate(vec2(pp2[0]+sideWidth*1.5, pp2[1]+sideWidth*0.35)), halign = RightAlign)
        image.fillText(font, labels[2], translate(vec2(pp3[0]-sideWidth*1.5, pp3[1]+sideWidth*0.35)), halign = LeftAlign)
    else:
        # If both bottom ones are chemical elements, we want them aligned along the axis
        image.fillText(font, labels[1], translate(vec2(pp2[0]+sideWidth*0.2, pp2[1]-sideWidth*0.3)), halign = LeftAlign)
        image.fillText(font, labels[2], translate(vec2(pp3[0]-sideWidth*0.2, pp3[1]-sideWidth*0.3)), halign = RightAlign)


# ********* Axis Ticks and Markers *********

proc drawAxisTicksMarkers(image: Image) =
    let ctx = newContext(image)
    ctx.fillStyle = rgba(5, 100, 80, 220)
    ctx.strokeStyle = rgba(5, 100, 80, 220)
    ctx.font = fontSupport
    ctx.fontsize = sideWidth*0.75

    if axisType == "decimal":
        # 10% increments
        ctx.textAlign = CenterAlign
        for i in 1..<10:
            ctx.fillText(
                "^", 
                vec2((pp2[0]*i.float + pp3[0]*(10-i.float))/10, pp2[1]+sideWidth*0.5+3*scaling))
        if longLabel:
            for i in 3..<8:
                ctx.fillText(
                    (10*i).intToStr, 
                    vec2((pp2[0]*i.float + pp3[0]*(10-i.float))/10, pp2[1]+sideWidth*1))
        else:
            for i in 1..<10:
                ctx.fillText(
                    (10*i).intToStr, 
                    vec2((pp2[0]*i.float + pp3[0]*(10-i.float))/10, pp2[1]+sideWidth*1))

        ctx.textAlign = LeftAlign
        for i in 1..<10:
            ctx.fillText(
                "<" & (10*i).intToStr, 
                vec2(
                    (pp1[0]*i.float + pp2[0]*(10-i.float))/10-sideWidth*0.050, 
                    (pp1[1]*i.float + pp2[1]*(10-i.float))/10+sideWidth*0.24)
            )

        ctx.textAlign = RightAlign
        for i in 1..<10:
            ctx.fillText(
                (10*(10-i)).intToStr & ">",
                vec2(
                    (pp1[0]*i.float + pp3[0]*(10-i.float))/10+sideWidth*0.050, 
                    (pp1[1]*i.float + pp3[1]*(10-i.float))/10+sideWidth*0.24)
            )
    elif axisType == "vulgarFractions":
        # 15 fractions like 1/8 or 2/3
        const 
            vulgarFractions: seq[string] = @["⅛", "⅙", "⅕", "¼", "⅓", "⅜", "⅖", "½", "⅗", "⅝", "⅔", "¾", "⅘", "⅚", "⅞"]
            vlugarPositions: seq[float] = @[1/8, 1/6, 1/5, 1/4, 1/3, 3/8, 2/5, 1/2, 3/5, 5/8, 2/3, 3/4, 4/5, 5/6, 7/8]
        
        ctx.textAlign = CenterAlign
        for i in 0..<15:
            let vp = vlugarPositions[i]
            ctx.fillText(
                "^", 
                vec2((pp2[0]*vp + pp3[0]*(1-vp)), pp2[1]+sideWidth*0.5+3*scaling))
        if longLabel:
            for i in 4..<11:
                let vp = vlugarPositions[i]
                ctx.fillText(
                    vulgarFractions[i], 
                    vec2((pp2[0]*vp + pp3[0]*(1-vp)), pp2[1]+sideWidth*1))
        else:
            for i in 0..<15:
                let vp = vlugarPositions[i]
                ctx.fillText(
                    vulgarFractions[i], 
                    vec2((pp2[0]*vp + pp3[0]*(1-vp)), pp2[1]+sideWidth*1))
        
        ctx.textAlign = LeftAlign
        for i in 0..<15:
            let vp = vlugarPositions[i]
            ctx.fillText(
                "<" & vulgarFractions[i], 
                vec2(
                    (pp1[0]*vp + pp2[0]*(1-vp))-sideWidth*0.050, 
                    (pp1[1]*vp + pp2[1]*(1-vp))+sideWidth*0.24)
            )
        
        ctx.textAlign = RightAlign
        for i in 0..<15:
            let vp = vlugarPositions[i]
            ctx.fillText(
                vulgarFractions[14-i] & ">",
                vec2(
                    (pp1[0]*vp + pp3[0]*(1-vp))+sideWidth*0.050, 
                    (pp1[1]*vp + pp3[1]*(1-vp))+sideWidth*0.24)
            )
    elif axisType == "quanta":
        # Marker at each step
        if ndiv<=24:
            ctx.textAlign = CenterAlign
            for i in 1..<ndiv:
                ctx.fillText(
                    "^", 
                    vec2((pp2[0]*i.float + pp3[0]*(ndiv-i).float)/ndiv.float, pp2[1]+sideWidth*0.5+3*scaling))
            if longLabel:
                if ndiv <= 12:
                    for i in 3..<ndiv-2:
                        ctx.fillText(
                            $i,
                            vec2((pp2[0]*i.float + pp3[0]*(ndiv-i).float)/ndiv.float, pp2[1]+sideWidth*1))
                else:
                    for i in 6..<ndiv-5:
                        ctx.fillText(
                            $i,
                            vec2((pp2[0]*i.float + pp3[0]*(ndiv-i).float)/ndiv.float, pp2[1]+sideWidth*1))
            else:
                for i in 1..<ndiv:
                    ctx.fillText(
                        $i, 
                        vec2((pp2[0]*i.float + pp3[0]*(ndiv-i).float)/ndiv.float, pp2[1]+sideWidth*1))
            ctx.textAlign = LeftAlign
            for i in 1..<ndiv:
                ctx.fillText(
                    "<" & $i, 
                    vec2(
                        (pp1[0]*i.float + pp2[0]*(ndiv-i).float)/ndiv.float-sideWidth*0.05,
                        (pp1[1]*i.float + pp2[1]*(ndiv-i).float)/ndiv.float+sideWidth*0.24)
                )
            ctx.textAlign = RightAlign
            for i in 1..<ndiv:
                ctx.fillText(
                    $i & ">", 
                    vec2(
                        (pp1[0]*i.float + pp3[0]*(ndiv-i).float)/ndiv.float+sideWidth*0.05,
                        (pp1[1]*i.float + pp3[1]*(ndiv-i).float)/ndiv.float+sideWidth*0.24)
                )

        else:
            ctx.fontsize = sideWidth*0.33
            ctx.textAlign = CenterAlign
            for i in 1..<ndiv:
                ctx.fillText(
                    "^", 
                    vec2((pp2[0]*i.float + pp3[0]*(ndiv-i).float)/ndiv.float, pp2[1]+sideWidth*0.2+3*scaling))
                ctx.fillText($i,vec2((pp2[0]*i.float + pp3[0]*(ndiv-i).float)/ndiv.float, pp2[1]+sideWidth*0.4))
            ctx.textAlign = LeftAlign
            for i in 1..<ndiv:
                ctx.fillText(
                    "<" & $i, 
                    vec2(
                        (pp1[0]*i.float + pp2[0]*(ndiv-i).float)/ndiv.float-sideWidth*0.03, 
                        (pp1[1]*i.float + pp2[1]*(ndiv-i).float)/ndiv.float+sideWidth*0.11)
                )
            ctx.textAlign = RightAlign
            for i in 1..<ndiv:
                ctx.fillText(
                    $i & ">", 
                    vec2(
                        (pp1[0]*i.float + pp3[0]*(ndiv-i).float)/ndiv.float+sideWidth*0.03, 
                        (pp1[1]*i.float + pp3[1]*(ndiv-i).float)/ndiv.float+sideWidth*0.11)
                )
            
                    


# ********* Elemental Legend *********
proc drawElementalLegend(image: Image) =
    let ctx = newContext(image)
    ctx.strokeStyle = rgba(0, 100, 100, 220)
    ctx.font = fontMain
    ctx.fontsize = sideWidth*1.5
    ctx.textAlign = LeftAlign

    const legendHexSize: float = (990/(nDiv+1)-1).clamp(15, 50) * scaling

    for (i, el) in elementList.pairs:
        let rgbMap = oklabToRGB((colorSchemeOKlab[i][0].float, colorSchemeOKlab[i][1].float, colorSchemeOKlab[i][2].float))

        ctx.fillStyle = rgba(rgbMap[0].uint8, rgbMap[1].uint8, rgbMap[2].uint8, alpha)
        ctx.fillPolygon(vec2(350*scaling, (170+100*i).float*scaling), legendHexSize, sides = 6)
        ctx.fillStyle = rgba(rgbMap[0].uint8, rgbMap[1].uint8, rgbMap[2].uint8, 255.uint8)
        ctx.fillText(
            el, 
            vec2(400*scaling, (100*i+200).float*scaling),
        )


# ********* Marker Legend *********
proc drawMarkerLegend(image: Image) = 
    let ctx = newContext(image)
    ctx.strokeStyle = rgba(0, 100, 100, 220)
    ctx.font = fontMain
    ctx.fontsize = sideWidth*1.5
    ctx.textAlign = LeftAlign
    if markerOverlay1:
        let position: float = 6.5
        ctx.fillStyle = rgba(0, 0, 0, 255)
        ctx.fillPolygon(vec2(350*scaling, (170+100*position).float*scaling), distance*1.05/2, sides = 24)
        ctx.fillStyle = rgba(200, 200, 0, 255)
        ctx.fillPolygon(vec2(350*scaling, (170+100*position).float*scaling), distance*1/2, sides = 24)
        ctx.fontsize = sideWidth
        ctx.fillText(
            marker1,
            vec2(400*scaling, (100*position+190).float*scaling),
        )

    if markerOverlay2:
        let position: float = 7.5
        ctx.fillStyle = rgba(0, 0, 0, 255)
        ctx.fillPolygon(vec2(350*scaling, (170+100*position).float*scaling), distance*0.75/2, sides = 24)
        ctx.fillStyle = rgba(255, 0, 0, 255)
        ctx.fillPolygon(vec2(350*scaling, (170+100*position).float*scaling), distance*0.7/2, sides = 24)
        ctx.fontsize = sideWidth
        ctx.fillText(
            marker2,
            vec2(400*scaling, (100*position+190).float*scaling),
        )


# ********* Property Legend *********
proc drawPropertyLegend(image: Image) =
    let ctx = newContext(image)
    ctx.strokeStyle = rgba(0, 100, 100, 220)
    ctx.font = fontMain
    ctx.fontsize = sideWidth
    ctx.textAlign = RightAlign

    # Property legend consisting of 100 levels (rectangles) and 11 line markers for 10% increments. Label at min and max.
    for i in 0..100:
        ctx.fillStyle = rgba(prop2rgb((100-i).float/100, propertyColoringOKlab))
        ctx.fillRect(2100*scaling, (120+5*i).float*scaling, 100*scaling, 5*scaling)
        if i mod 10 == 0:
            ctx.fillRect(2050*scaling, (120+5.03*i.float)*scaling, 125*scaling, 2*scaling)

    # Max label
    ctx.fillStyle = rgba(0, 0, 0, 255)
    for dx in -1..1:
        for dy in -1..1:
            ctx.fillText(
                $maxVal, 
                vec2((2020.float+dx.float)*scaling, (120.float+dy.float)*scaling+sideWidth/2),
            )
    ctx.fillStyle = prop2rgb(1, propertyColoringOKlab)
    ctx.fillText(
        $maxVal, 
        vec2(2020*scaling, 120*scaling+sideWidth/2),
    )

    # Min label
    ctx.fillStyle = rgba(0, 0, 0, 255)
    for dx in -1..1:
        for dy in -1..1:
            ctx.fillText(
                $minVal, 
                vec2((2020.float+dx.float)*scaling, (620.float+dy.float)*scaling+sideWidth/2),
            )
    ctx.fillStyle = prop2rgb(0, propertyColoringOKlab)
    ctx.fillText(
        $minVal, 
        vec2(2020*scaling, 620*scaling+sideWidth/2),
    )

# ********* Index Overlay *********
proc drawIndicies(image: Image): void =
    let ctx = newContext(image)
    ctx.fillStyle = rgba(75, 75, 75, 255)
    ctx.font = fontSupport
    ctx.fontsize = thinLine*4
    ctx.textAlign = CenterAlign
    for i in 0..<gpl.len:
        let p = gpl[i]
        ctx.fillText($i, vec2(p[0], p[1]-thinLine*2))

# ********* Final Processing *********

let image = newImage(width, height)
image.fillWhite()
image.drawBackground()
image.drawCompositonHexes(elementalColoring)
if propertyOverlay: 
    image.drawPropertyHexes(propertyColoring)
if pathPointsOverlay:
    image.drawDesignedPath(pathPoints)
if markerOverlay1: 
    image.drawMarkers1(feasibilityField1)
if markerOverlay2: 
    image.drawMarkers2(feasibilityField1, feasibilityField2)
image.drawForeground()
image.drawAxisLabels()
image.drawAxisTicksMarkers()
image.drawElementalLegend()
if markerOverlay1 or markerOverlay2:
    image.drawMarkerLegend()
if propertyOverlay:
    image.drawPropertyLegend()
if indexOverlay:
    image.drawIndicies()

# ********* Save the image *********
let t1 = cpuTime()
image.writeFile(filename)

let nPrimitives: int = (gpl.len*(2+2) + 10 + 30)
echo "Approximately ", nPrimitives, " primitives"
echo "Constructed in ", (t1 - t0).round(3), "s and processed into PNG in ", (cpuTime() - t1).round(3), "s"