var Vector = function(x, y) {
	this.x = x;
	this.y = y;
};

Vector.polar = function(r, theta) {
	return new Vector(r * Math.cos(theta), r * Math.sin(theta));
};

Vector.prototype.length = function() {
	return Math.sqrt(this.x * this.x + this.y * this.y);
};

Vector.prototype.copy = function() {
	return new Vector(this.x, this.y);
};

Vector.prototype.setLength = function(length) {
	var current = this.length();
	if (current) {
		var scale = length / current;
		this.x *= scale;
		this.y *= scale;
	}
	return this;
};

Vector.prototype.setAngle = function(theta) {
	var r = length();
	this.x = r * Math.cos(theta);
	this.y = r * Math.sin(theta);
	return this;
};

Vector.prototype.getAngle = function() {
	return Math.atan2(this.y, this.x);
};

Vector.prototype.d = function(v) {
	var dx = v.x - this.x;
	var dy = v.y - this.y;
	return Math.sqrt(dx * dx + dy * dy);
};

var IDProjection = {

	project: function(x, y, opt_v) {
		var v = opt_v || new Vector();
		v.x = x;
		v.y = y;
		return v;
	},

	invert: function(x, y, opt_v) {
		var v = opt_v || new Vector();
		v.x = x;
		v.y = y;
		return v;
	}
};

var Albers = function() {

	function radians(degrees) {
		return Math.PI * degrees / 180;
	}

	var phi1 = radians(29.5);
	var phi2 = radians(45.5);
	var n = 0.5 * (phi1 + phi2);
	var C = Math.cos(phi1) * Math.cos(phi1) + 2 * n * Math.sin(phi1);
	var phi0 = radians(38);
	var lambda0 = radians(-98);
	var rho0 = Math.sqrt(C - 2 * n * Math.sin(phi0)) / n;

	return {

		project: function(lon, lat, opt_result) {
			lon = radians(lon);
			lat = radians(lat);
			var theta = n * (lon - lambda0);
			var rho = Math.sqrt(C - 2 * n * Math.sin(lat)) / n;
			var x = rho * Math.sin(theta);
			var y = rho0 - rho * Math.cos(theta);
			if (opt_result) {
				opt_result.x = x;
				opt_result.y = y;
				return opt_result;
			}
			return new Vector(x, y);
		},

		invert: function(x, y) {
			var rho2 = x * x + (rho0 - y) * (rho0 - y);
			var theta = Math.atan(x / (rho0 - y));
			var lon = lambda0 + theta / n;
			var lat = Math.asin((C / n - rho2 * n) / 2);
			return new Vector(lon * 180 / Math.PI, lat * 180 / Math.PI);
		}
	};
}();

var ScaledAlbers = function(scale, offsetX, offsetY, longMin, latMin) {
	this.scale = scale;
	this.offsetX = offsetX;
	this.offsetY = offsetY;
	this.longMin = longMin;
	this.latMin = latMin;
	this.swCorner = Albers.project(longMin, latMin);
};

ScaledAlbers.temp = new Vector(0, 0);

ScaledAlbers.prototype.project = function(lon, lat, opt_result) {
	var proj = Albers.project(lon, lat, ScaledAlbers.temp);
	var a = proj.x;
	var b = proj.y;
	var x = this.scale * (a - this.swCorner.x) + this.offsetX;
	var y = -this.scale * (b - this.swCorner.y) + this.offsetY;
	if (opt_result) {
		opt_result.x = x;
		opt_result.y = y;
		return opt_result;
	}
	return new Vector(x, y);
};

ScaledAlbers.prototype.invert = function(x, y) {
	var a = (x - this.offsetX) / this.scale + this.swCorner.x;
	var b = (y - this.offsetY) / -this.scale + this.swCorner.y;
	return Albers.invert(a, b);
};

var VectorField = function(field, x0, y0, x1, y1) {
	this.x0 = x0;
	this.x1 = x1;
	this.y0 = y0;
	this.y1 = y1;
	this.field = field;
	this.w = field.length;
	this.h = field[0].length;
	this.maxLength = 0;
	var mx = 0;
	var my = 0;
	for (var i = 0; i < this.w; i++) {
		for (var j = 0; j < this.h; j++) {
			if (field[i][j].length() > this.maxLength) {
				mx = i;
				my = j;
			}
			this.maxLength = Math.max(this.maxLength, field[i][j].length());
		}
	}
	mx = (mx / this.w) * (x1 - x0) + x0;
	my = (my / this.h) * (y1 - y0) + y0;
};

VectorField.read = function(data, correctForSphere) {
	var field = [];
	var w = data.gridWidth;
	var h = data.gridHeight;
	var n = 2 * w * h;
	var i = 0;

	var total = 0;
	var weight = 0;

	for (var x = 0; x < w; x++) {
		field[x] = [];
		for (var y = 0; y < h; y++) {
			var vx = data.field[i++];
			var vy = data.field[i++];
			var v = new Vector(vx, vy);
			// Test vector below:
			// v = new Vector(10, -10);
			if (correctForSphere) {
				var ux = x / (w - 1);
				var uy = y / (h - 1);
				var lon = data.x0 * (1 - ux) + data.x1 * ux;
				var lat = data.y0 * (1 - uy) + data.y1 * uy;
				var m = Math.PI * lat / 180;
				var length = v.length();
				if (length) {
					total += length * m;
					weight += m;
				}
				v.x /= Math.cos(m);
				v.setLength(length);
			}
			field[x][y] = v;
		}
	}
	var result = new VectorField(field, data.x0, data.y0, data.x1, data.y1);

	if (total && weight) {
		result.averageLength = total / weight;
	}
	return result;
};

VectorField.prototype.inBounds = function(x, y) {
	return x >= this.x0 && x < this.x1 && y >= this.y0 && y < this.y1;
};

VectorField.prototype.bilinear = function(coord, a, b) {
	var na = Math.floor(a);
	var nb = Math.floor(b);
	var ma = Math.ceil(a);
	var mb = Math.ceil(b);
	var fa = a - na;
	var fb = b - nb;

	return	this.field[na][nb][coord] * (1 - fa) * (1 - fb) +
			this.field[ma][nb][coord] * fa * (1 - fb) +
			this.field[na][mb][coord] * (1 - fa) * fb +
			this.field[ma][mb][coord] * fa * fb;
};

VectorField.prototype.getValue = function(x, y, opt_result) {
	var a = (this.w - 1 - 1e-6) * (x - this.x0) / (this.x1 - this.x0);
	var b = (this.h - 1 - 1e-6) * (y - this.y0) / (this.y1 - this.y0);
	var vx = this.bilinear('x', a, b);
	var vy = this.bilinear('y', a, b);
	if (opt_result) {
		opt_result.x = vx;
		opt_result.y = vy;
		return opt_result;
	}
	return new Vector(vx, vy);
};

VectorField.prototype.vectValue = function(vector) {
	return this.getValue(vector.x, vector.y);
};

VectorField.constant = function(dx, dy, x0, y0, x1, y1) {
	var field = new VectorField([[]], x0, y0, x1, y1);
	field.maxLength = Math.sqrt(dx * dx + dy * dy);
	field.getValue = function() {
		return new Vector(dx, dy);
	}
	return field;
};

var Animator = function(element, opt_animFunc, opt_unzoomButton) {
	this.element = element;
	this.mouseIsDown = false;
	this.mouseX = -1;
	this.mouseY = -1;
	this.animating = true;
	this.state = 'animate';
	this.listeners = [];
	this.dx = 0;
	this.dy = 0;
	this.scale = 1;
	this.zoomProgress = 0;
	this.scaleTarget = 1;
	this.scaleStart = 1;
	this.animFunc = opt_animFunc;
	this.unzoomButton = opt_unzoomButton;

	if (element) {
		var self = this;
		$(element).mousedown(function(e) {
			self.mouseX = e.pageX - this.offsetLeft;
			self.mouseY = e.pageY - this.offsetTop;
			self.mousedown();
		});
		$(element).mouseup(function(e) {
			self.mouseX = e.pageX - this.offsetLeft;
			self.mouseY = e.pageY - this.offsetTop;
			self.mouseup();
		});
		$(element).mousemove(function(e) {
			self.mouseX = e.pageX - this.offsetLeft;
			self.mouseY = e.pageY - this.offsetTop;
			self.mousemove();
		});
	}
};

Animator.prototype.mousedown = function() {
	this.state = 'mouse-down';
	this.notify('startMove');
	this.landingX = this.mouseX;
	this.landingY = this.mouseY;
	this.dxStart = this.dx;
	this.dyStart = this.dy;
	this.scaleStart = this.scale;
	this.mouseIsDown = true;
};

Animator.prototype.mousemove = function() {
	if (!this.mouseIsDown) {
		this.notify('hover');
		return;
	}
	var ddx = this.mouseX - this.landingX;
	var ddy = this.mouseY - this.landingY;
	var slip = Math.abs(ddx) + Math.abs(ddy);
	if (slip > 2 || this.state == 'pan') {
		this.state = 'pan';
		this.dx += ddx;
		this.dy += ddy;
		this.landingX = this.mouseX;
		this.landingY = this.mouseY;
		this.notify('move');
	}
};

Animator.prototype.mouseup = function() {
	this.mouseIsDown = false;
	if (this.state == 'pan') {
		this.state = 'animate';
		this.notify('endMove');
		return;
	}
	this.zoomClick(this.mouseX, this.mouseY);
};

Animator.prototype.add = function(listener) {
	this.listeners.push(listener);
};

Animator.prototype.notify = function(message) {
	if (this.unzoomButton) {
		var diff = Math.abs(this.scale - 1) > 0.001 || Math.abs(this.dx) > 0.001 || Math.abs(this.dy > 0.001);
		this.unzoomButton.style.visibility = diff ? 'visible' : 'hidden';
	}
	if (this.animFunc && !this.animFunc()) {
		return;
	}
	for (var i = 0; i < this.listeners.length; i++) {
		var listener = this.listeners[i];
		if (listener[message]) {
			listener[message].call(listener, this);
		}
	}
};

Animator.prototype.unzoom = function() {
	this.zoom(0, 0, 1);
};

Animator.prototype.zoomClick = function(x, y) {
	var z = 1.7;
	var scale = 1.7 * this.scale;
	var dx = x - z * (x - this.dx);
	var dy = y - z * (y - this.dy);
	this.zoom(dx, dy, scale);
};

Animator.prototype.zoom = function(dx, dy, scale) {
	this.state = 'zoom';
	this.zoomProgress = 0;
	this.scaleStart = this.scale;
	this.scaleTarget = scale;
	this.dxTarget = dx;
	this.dyTarget = dy;
	this.dxStart = this.dx;
	this.dyStart = this.dy;
	this.notify('startMove');
};

Animator.prototype.relativeZoom = function() {
	return this.scale / this.scaleStart;
};

Animator.prototype.relativeDx = function() {
	return this.dx - this.dxStart;
};

Animator.prototype.relativeDy = function() {
	return this.dy - this.dyStart;
};

Animator.prototype.start = function(opt_millis) {
	var millis = opt_millis || 20;
	var self = this;

	function go() {
		var start = new Date();
		self.loop();
		var time = new Date() - start;
		setTimeout(go, Math.max(10, millis - time));
	}
	go();
};

Animator.prototype.loop = function() {
	if (this.state == 'mouse-down' || this.state == 'pan') {
		return;
	}

	if (this.state == 'animate') {
		this.notify('animate');
		return;
	}

	if (this.state == 'zoom') {
		this.zoomProgress = Math.min(1, this.zoomProgress + 0.07);
		var u = (1 + Math.cos(Math.PI * this.zoomProgress)) / 2;
		
		function lerp(a, b) {
			return u * a + (1 - u) * b;
		}
		this.scale = lerp(this.scaleStart, this.scaleTarget);
		this.dx = lerp(this.dxStart, this.dxTarget);
		this.dy = lerp(this.dyStart, this.dyTarget);

		if (this.zoomProgress < 1) {
			this.notify('move');
		}
		else {
			this.state = 'animate';
			this.zoomCurrent = this.zoomTarget;
			this.notify('endMove');
		}
	}
};

var Particle = function(x, y, age) {
	this.x = x;
	this.y = y;
	this.oldX = -1;
	this.oldY = -1;
	this.age = age;
	this.rnd = Math.random();
};

var MotionDisplay = function(canvas, imageCanvas, field, numParticles, opt_projection) {
	this.canvas = canvas;
	this.projection = opt_projection || IDProjection;
	this.field = field;
	this.numParticles = numParticles;
	this.first = true;
	this.maxLength = field.maxLength;
	this.speedScale = 0.5;
	this.renderState = 'normal';
	this.imageCanvas = imageCanvas;
	this.x0 = this.field.x0;
	this.x1 = this.field.x1;
	this.y0 = this.field.y0;
	this.y1 = this.field.y1;
	this.makeNewParticles(null, true);
	this.colors = [];
	this.rgb = '40, 40, 40';
	this.background = 'rgb(' + this.rgb + ')';
	this.backgroundAlpha = 'rgba(' + this.rgb + ', 0.02)';
	this.outsideColor = '#FFFFFF';
	// this.colors[0] = 'rgb(111, 112, 114)';
	// this.colors[1] = 'rgb(111, 112, 114)';
	// this.colors[2] = 'rgb(105, 119, 121)';
	// this.colors[3] = 'rgb(98, 127, 128)';
	// this.colors[4] = 'rgb(92, 134, 135)';
	// this.colors[5] = 'rgb(85, 141, 142)';
	// this.colors[6] = 'rgb(79, 149, 149)';
	// this.colors[7] = 'rgb(73, 156, 156)';
	// this.colors[8] = 'rgb(66, 163, 163)';
	// this.colors[9] = 'rgb(60, 171, 170)';
	// this.colors[10] = 'rgb(53, 178, 177)';
	// this.colors[11] = 'rgb(47, 186, 184)';
	// this.colors[12] = 'rgb(40, 193, 191)';
	// this.colors[13] = 'rgb(34, 200, 198)';
	// this.colors[14] = 'rgb(28, 208, 205)';
	// this.colors[15] = 'rgb(21, 215, 212)';
	// this.colors[16] = 'rgb(15, 222, 219)';
	// this.colors[17] = 'rgb(8, 230, 226)';
	// this.colors[18] = 'rgb(2, 237, 233)';
	// this.colors[19] = 'rgb(2, 233, 234)';
	// this.colors[20] = 'rgb(2, 229, 235)';
	// this.colors[21] = 'rgb(2, 225, 236)';
	// this.colors[22] = 'rgb(2, 221, 237)';
	// this.colors[23] = 'rgb(3, 217, 237)';
	// this.colors[24] = 'rgb(3, 213, 238)';
	// this.colors[25] = 'rgb(3, 209, 239)';
	// this.colors[26] = 'rgb(3, 205, 240)';
	// this.colors[27] = 'rgb(3, 200, 241)';
	// this.colors[28] = 'rgb(3, 196, 242)';
	// this.colors[29] = 'rgb(3, 192, 243)';
	// this.colors[30] = 'rgb(3, 188, 244)';
	// this.colors[31] = 'rgb(4, 184, 244)';
	// this.colors[32] = 'rgb(4, 180, 245)';
	// this.colors[33] = 'rgb(4, 176, 246)';
	// this.colors[34] = 'rgb(4, 172, 247)';
	// this.colors[35] = 'rgb(4, 168, 248)';
	// this.colors[36] = 'rgb(5, 158, 248)';
	// this.colors[37] = 'rgb(5, 148, 248)';
	// this.colors[38] = 'rgb(6, 138, 249)';
	// this.colors[39] = 'rgb(7, 128, 249)';
	// this.colors[40] = 'rgb(7, 119, 249)';
	// this.colors[41] = 'rgb(8, 109, 249)';
	// this.colors[42] = 'rgb(9, 99, 250)';
	// this.colors[43] = 'rgb(9, 89, 250)';
	// this.colors[44] = 'rgb(10, 79, 250)';
	// this.colors[45] = 'rgb(10, 69, 250)';
	// this.colors[46] = 'rgb(11, 59, 251)';
	// this.colors[47] = 'rgb(12, 49, 251)';
	// this.colors[48] = 'rgb(12, 40, 251)';
	// this.colors[49] = 'rgb(13, 30, 251)';
	// this.colors[50] = 'rgb(14, 20, 252)';
	// this.colors[51] = 'rgb(14, 10, 252)';
	// this.colors[52] = 'rgb(15, 0, 252)';
	// this.colors[53] = 'rgb(14, 15, 237)';
	// this.colors[54] = 'rgb(13, 30, 222)';
	// this.colors[55] = 'rgb(12, 45, 208)';
	// this.colors[56] = 'rgb(11, 60, 193)';
	// this.colors[57] = 'rgb(11, 75, 178)';
	// this.colors[58] = 'rgb(10, 90, 163)';
	// this.colors[59] = 'rgb(9, 105, 148)';
	// this.colors[60] = 'rgb(8, 120, 133)';
	// this.colors[61] = 'rgb(7, 134, 119)';
	// this.colors[62] = 'rgb(6, 149, 104)';
	// this.colors[63] = 'rgb(5, 164, 89)';
	// this.colors[64] = 'rgb(4, 179, 74)';
	// this.colors[65] = 'rgb(4, 194, 59)';
	// this.colors[66] = 'rgb(3, 209, 44)';
	// this.colors[67] = 'rgb(2, 224, 30)';
	// this.colors[68] = 'rgb(1, 239, 15)';
	// this.colors[69] = 'rgb(0, 254, 0)';
	// this.colors[70] = 'rgb(0, 251, 0)';
	// this.colors[71] = 'rgb(0, 248, 0)';
	// this.colors[72] = 'rgb(0, 246, 0)';
	// this.colors[73] = 'rgb(0, 243, 0)';
	// this.colors[74] = 'rgb(0, 240, 0)';
	// this.colors[75] = 'rgb(0, 237, 0)';
	// this.colors[76] = 'rgb(0, 234, 0)';
	// this.colors[77] = 'rgb(0, 231, 0)';
	// this.colors[78] = 'rgb(0, 229, 0)';
	// this.colors[79] = 'rgb(0, 226, 0)';
	// this.colors[80] = 'rgb(0, 223, 0)';
	// this.colors[81] = 'rgb(0, 220, 0)';
	// this.colors[82] = 'rgb(0, 217, 0)';
	// this.colors[83] = 'rgb(0, 214, 0)';
	// this.colors[84] = 'rgb(0, 212, 0)';
	// this.colors[85] = 'rgb(0, 209, 0)';
	// this.colors[86] = 'rgb(0, 206, 0)';
	// this.colors[87] = 'rgb(0, 203, 0)';
	// this.colors[88] = 'rgb(0, 200, 0)';
	// this.colors[89] = 'rgb(0, 197, 0)';
	// this.colors[90] = 'rgb(0, 194, 0)';
	// this.colors[91] = 'rgb(0, 190, 0)';
	// this.colors[92] = 'rgb(0, 187, 0)';
	// this.colors[93] = 'rgb(0, 184, 0)';
	// this.colors[94] = 'rgb(0, 181, 0)';
	// this.colors[95] = 'rgb(0, 178, 0)';
	// this.colors[96] = 'rgb(0, 175, 0)';
	// this.colors[97] = 'rgb(0, 172, 0)';
	// this.colors[98] = 'rgb(0, 169, 0)';
	// this.colors[99] = 'rgb(0, 165, 0)';
	// this.colors[100] = 'rgb(0, 162, 0)';
	// this.colors[101] = 'rgb(0, 159, 0)';
	// this.colors[102] = 'rgb(0, 156, 0)';
	// this.colors[103] = 'rgb(0, 153, 0)';
	// this.colors[104] = 'rgb(15, 159, 0)';
	// this.colors[105] = 'rgb(30, 164, 0)';
	// this.colors[106] = 'rgb(45, 170, 0)';
	// this.colors[107] = 'rgb(60, 175, 0)';
	// this.colors[108] = 'rgb(75, 181, 0)';
	// this.colors[109] = 'rgb(90, 187, 0)';
	// this.colors[110] = 'rgb(105, 192, 0)';
	// this.colors[111] = 'rgb(120, 198, 0)';
	// this.colors[112] = 'rgb(135, 203, 0)';
	// this.colors[113] = 'rgb(150, 209, 0)';
	// this.colors[114] = 'rgb(165, 214, 0)';
	// this.colors[115] = 'rgb(180, 220, 0)';
	// this.colors[116] = 'rgb(195, 226, 0)';
	// this.colors[117] = 'rgb(210, 231, 0)';
	// this.colors[118] = 'rgb(225, 237, 0)';
	// this.colors[119] = 'rgb(240, 242, 0)';
	// this.colors[120] = 'rgb(255, 248, 0)';
	// this.colors[121] = 'rgb(254, 245, 0)';
	// this.colors[122] = 'rgb(252, 242, 0)';
	// this.colors[123] = 'rgb(251, 239, 0)';
	// this.colors[124] = 'rgb(249, 236, 0)';
	// this.colors[125] = 'rgb(248, 232, 0)';
	// this.colors[126] = 'rgb(247, 229, 0)';
	// this.colors[127] = 'rgb(245, 226, 0)';
	// this.colors[128] = 'rgb(244, 223, 0)';
	// this.colors[129] = 'rgb(242, 220, 0)';
	// this.colors[130] = 'rgb(241, 217, 0)';
	// this.colors[131] = 'rgb(239, 214, 0)';
	// this.colors[132] = 'rgb(238, 211, 0)';
	// this.colors[133] = 'rgb(237, 207, 0)';
	// this.colors[134] = 'rgb(235, 204, 0)';
	// this.colors[135] = 'rgb(234, 201, 0)';
	// this.colors[136] = 'rgb(232, 198, 0)';
	// this.colors[137] = 'rgb(231, 195, 0)';
	// this.colors[138] = 'rgb(232, 193, 0)';
	// this.colors[139] = 'rgb(233, 192, 1)';
	// this.colors[140] = 'rgb(235, 190, 1)';
	// this.colors[141] = 'rgb(236, 189, 2)';
	// this.colors[142] = 'rgb(237, 187, 2)';
	// this.colors[143] = 'rgb(238, 185, 3)';
	// this.colors[144] = 'rgb(239, 184, 3)';
	// this.colors[145] = 'rgb(240, 182, 4)';
	// this.colors[146] = 'rgb(242, 181, 4)';
	// this.colors[147] = 'rgb(243, 179, 5)';
	// this.colors[148] = 'rgb(244, 178, 5)';
	// this.colors[149] = 'rgb(245, 176, 6)';
	// this.colors[150] = 'rgb(246, 174, 6)';
	// this.colors[151] = 'rgb(247, 173, 7)';
	// this.colors[152] = 'rgb(249, 171, 7)';
	// this.colors[153] = 'rgb(250, 170, 8)';
	// this.colors[154] = 'rgb(251, 168, 8)';
	// this.colors[155] = 'rgb(251, 158, 8)';
	// this.colors[156] = 'rgb(251, 149, 7)';
	// this.colors[157] = 'rgb(251, 139, 7)';
	// this.colors[158] = 'rgb(251, 129, 6)';
	// this.colors[159] = 'rgb(252, 119, 6)';
	// this.colors[160] = 'rgb(252, 110, 5)';
	// this.colors[161] = 'rgb(252, 100, 5)';
	// this.colors[162] = 'rgb(252, 90, 4)';
	// this.colors[163] = 'rgb(252, 81, 4)';
	// this.colors[164] = 'rgb(252, 71, 3)';
	// this.colors[165] = 'rgb(252, 61, 3)';
	// this.colors[166] = 'rgb(252, 52, 2)';
	// this.colors[167] = 'rgb(253, 42, 2)';
	// this.colors[168] = 'rgb(253, 32, 1)';
	// this.colors[169] = 'rgb(253, 22, 1)';
	// this.colors[170] = 'rgb(253, 13, 0)';
	// this.colors[171] = 'rgb(253, 3, 0)';
	// this.colors[172] = 'rgb(251, 3, 0)';
	// this.colors[173] = 'rgb(250, 3, 0)';
	// this.colors[174] = 'rgb(248, 3, 0)';
	// this.colors[175] = 'rgb(246, 3, 0)';
	// this.colors[176] = 'rgb(244, 3, 0)';
	// this.colors[177] = 'rgb(243, 3, 0)';
	// this.colors[178] = 'rgb(241, 3, 0)';
	// this.colors[179] = 'rgb(239, 3, 0)';
	// this.colors[180] = 'rgb(238, 3, 0)';
	// this.colors[181] = 'rgb(236, 3, 0)';
	// this.colors[182] = 'rgb(234, 3, 0)';
	// this.colors[183] = 'rgb(233, 3, 0)';
	// this.colors[184] = 'rgb(231, 3, 0)';
	// this.colors[185] = 'rgb(229, 3, 0)';
	// this.colors[186] = 'rgb(227, 3, 0)';
	// this.colors[187] = 'rgb(226, 3, 0)';
	// this.colors[188] = 'rgb(224, 3, 0)';
	// this.colors[189] = 'rgb(222, 3, 0)';
	// this.colors[190] = 'rgb(221, 3, 0)';
	// this.colors[191] = 'rgb(219, 3, 0)';
	// this.colors[192] = 'rgb(217, 3, 0)';
	// this.colors[193] = 'rgb(216, 3, 0)';
	// this.colors[194] = 'rgb(214, 3, 0)';
	// this.colors[195] = 'rgb(212, 3, 0)';
	// this.colors[196] = 'rgb(211, 3, 0)';
	// this.colors[197] = 'rgb(209, 2, 0)';
	// this.colors[198] = 'rgb(208, 2, 0)';
	// this.colors[199] = 'rgb(206, 2, 0)';
	// this.colors[200] = 'rgb(204, 2, 0)';
	// this.colors[201] = 'rgb(203, 2, 0)';
	// this.colors[202] = 'rgb(201, 2, 0)';
	// this.colors[203] = 'rgb(199, 2, 0)';
	// this.colors[204] = 'rgb(198, 2, 0)';
	// this.colors[205] = 'rgb(196, 2, 0)';
	// this.colors[206] = 'rgb(199, 2, 15)';
	// this.colors[207] = 'rgb(201, 2, 30)';
	// this.colors[208] = 'rgb(204, 2, 45)';
	// this.colors[209] = 'rgb(207, 2, 60)';
	// this.colors[210] = 'rgb(209, 1, 75)';
	// this.colors[211] = 'rgb(212, 1, 90)';
	// this.colors[212] = 'rgb(215, 1, 105)';
	// this.colors[213] = 'rgb(217, 1, 120)';
	// this.colors[214] = 'rgb(220, 1, 135)';
	// this.colors[215] = 'rgb(222, 1, 150)';
	// this.colors[216] = 'rgb(225, 1, 165)';
	// this.colors[217] = 'rgb(228, 1, 180)';
	// this.colors[218] = 'rgb(230, 0, 195)';
	// this.colors[219] = 'rgb(233, 0, 210)';
	// this.colors[220] = 'rgb(236, 0, 225)';
	// this.colors[221] = 'rgb(238, 0, 240)';
	// this.colors[222] = 'rgb(241, 0, 255)';
	// this.colors[223] = 'rgb(236, 5, 252)';
	// this.colors[224] = 'rgb(231, 11, 249)';
	// this.colors[225] = 'rgb(227, 16, 246)';
	// this.colors[226] = 'rgb(222, 22, 243)';
	// this.colors[227] = 'rgb(217, 27, 240)';
	// this.colors[228] = 'rgb(212, 33, 237)';
	// this.colors[229] = 'rgb(208, 38, 234)';
	// this.colors[230] = 'rgb(203, 44, 231)';
	// this.colors[231] = 'rgb(198, 49, 227)';
	// this.colors[232] = 'rgb(193, 55, 224)';
	// this.colors[233] = 'rgb(189, 60, 221)';
	// this.colors[234] = 'rgb(184, 66, 218)';
	// this.colors[235] = 'rgb(179, 71, 215)';
	// this.colors[236] = 'rgb(174, 77, 212)';
	// this.colors[237] = 'rgb(170, 82, 209)';
	// this.colors[238] = 'rgb(165, 88, 206)';
	// this.colors[239] = 'rgb(160, 93, 203)';
	// this.colors[240] = 'rgb(166, 102, 206)';
	// this.colors[241] = 'rgb(171, 112, 209)';
	// this.colors[242] = 'rgb(177, 121, 212)';
	// this.colors[243] = 'rgb(182, 131, 215)';
	// this.colors[244] = 'rgb(188, 140, 218)';
	// this.colors[245] = 'rgb(194, 150, 221)';
	// this.colors[246] = 'rgb(199, 159, 224)';
	// this.colors[247] = 'rgb(205, 169, 227)';
	// this.colors[248] = 'rgb(210, 178, 231)';
	// this.colors[249] = 'rgb(216, 188, 234)';
	// this.colors[250] = 'rgb(221, 197, 237)';
	// this.colors[251] = 'rgb(227, 207, 240)';
	// this.colors[252] = 'rgb(233, 216, 243)';
	// this.colors[253] = 'rgb(238, 226, 246)';
	// this.colors[254] = 'rgb(244, 235, 249)';
	// this.colors[255] = 'rgb(249, 245, 252)';
	this.colors[0] = 'rgb(232, 151, 254)';
	this.colors[1] = 'rgb(230, 143, 254)';
	this.colors[2] = 'rgb(228, 135, 253)';
	this.colors[3] = 'rgb(226, 126, 253)';
	this.colors[4] = 'rgb(225, 118, 253)';
	this.colors[5] = 'rgb(223, 109, 253)';
	this.colors[6] = 'rgb(221, 101, 253)';
	this.colors[7] = 'rgb(219, 93, 253)';
	this.colors[8] = 'rgb(216, 84, 252)';
	this.colors[9] = 'rgb(213, 76, 252)';
	this.colors[10] = 'rgb(210, 67, 252)';
	this.colors[11] = 'rgb(207, 59, 252)';
	this.colors[12] = 'rgb(204, 51, 251)';
	this.colors[13] = 'rgb(201, 42, 251)';
	this.colors[14] = 'rgb(197, 35, 251)';
	this.colors[15] = 'rgb(193, 30, 250)';
	this.colors[16] = 'rgb(188, 25, 250)';
	this.colors[17] = 'rgb(184, 20, 249)';
	this.colors[18] = 'rgb(179, 15, 248)';
	this.colors[19] = 'rgb(175, 10, 248)';
	this.colors[20] = 'rgb(170, 5, 247)';
	this.colors[21] = 'rgb(166, 1, 247)';
	this.colors[22] = 'rgb(162, 1, 246)';
	this.colors[23] = 'rgb(158, 1, 245)';
	this.colors[24] = 'rgb(154, 0, 244)';
	this.colors[25] = 'rgb(150, 0, 243)';
	this.colors[26] = 'rgb(145, 0, 242)';
	this.colors[27] = 'rgb(141, 0, 241)';
	this.colors[28] = 'rgb(138, 0, 240)';
	this.colors[29] = 'rgb(134, 0, 239)';
	this.colors[30] = 'rgb(131, 0, 238)';
	this.colors[31] = 'rgb(128, 0, 238)';
	this.colors[32] = 'rgb(124, 0, 237)';
	this.colors[33] = 'rgb(121, 0, 236)';
	this.colors[34] = 'rgb(118, 0, 235)';
	this.colors[35] = 'rgb(115, 1, 234)';
	this.colors[36] = 'rgb(113, 4, 234)';
	this.colors[37] = 'rgb(111, 7, 233)';
	this.colors[38] = 'rgb(108, 9, 232)';
	this.colors[39] = 'rgb(106, 12, 232)';
	this.colors[40] = 'rgb(104, 15, 231)';
	this.colors[41] = 'rgb(102, 18, 230)';
	this.colors[42] = 'rgb(100, 22, 230)';
	this.colors[43] = 'rgb(99, 26, 230)';
	this.colors[44] = 'rgb(98, 31, 231)';
	this.colors[45] = 'rgb(97, 35, 231)';
	this.colors[46] = 'rgb(95, 40, 231)';
	this.colors[47] = 'rgb(94, 44, 231)';
	this.colors[48] = 'rgb(93, 49, 231)';
	this.colors[49] = 'rgb(92, 54, 232)';
	this.colors[50] = 'rgb(91, 59, 232)';
	this.colors[51] = 'rgb(91, 64, 233)';
	this.colors[52] = 'rgb(90, 69, 234)';
	this.colors[53] = 'rgb(89, 74, 235)';
	this.colors[54] = 'rgb(89, 79, 236)';
	this.colors[55] = 'rgb(88, 84, 236)';
	this.colors[56] = 'rgb(87, 89, 238)';
	this.colors[57] = 'rgb(87, 95, 239)';
	this.colors[58] = 'rgb(86, 101, 240)';
	this.colors[59] = 'rgb(86, 107, 242)';
	this.colors[60] = 'rgb(85, 113, 243)';
	this.colors[61] = 'rgb(85, 118, 244)';
	this.colors[62] = 'rgb(84, 124, 246)';
	this.colors[63] = 'rgb(84, 131, 247)';
	this.colors[64] = 'rgb(84, 137, 248)';
	this.colors[65] = 'rgb(83, 143, 249)';
	this.colors[66] = 'rgb(83, 149, 251)';
	this.colors[67] = 'rgb(83, 156, 252)';
	this.colors[68] = 'rgb(82, 162, 253)';
	this.colors[69] = 'rgb(82, 168, 254)';
	this.colors[70] = 'rgb(81, 174, 254)';
	this.colors[71] = 'rgb(81, 179, 254)';
	this.colors[72] = 'rgb(80, 185, 254)';
	this.colors[73] = 'rgb(79, 190, 254)';
	this.colors[74] = 'rgb(79, 196, 254)';
	this.colors[75] = 'rgb(78, 202, 254)';
	this.colors[76] = 'rgb(78, 207, 254)';
	this.colors[77] = 'rgb(76, 211, 252)';
	this.colors[78] = 'rgb(75, 215, 249)';
	this.colors[79] = 'rgb(74, 219, 247)';
	this.colors[80] = 'rgb(73, 223, 245)';
	this.colors[81] = 'rgb(72, 227, 243)';
	this.colors[82] = 'rgb(71, 231, 241)';
	this.colors[83] = 'rgb(69, 235, 239)';
	this.colors[84] = 'rgb(67, 238, 235)';
	this.colors[85] = 'rgb(65, 240, 231)';
	this.colors[86] = 'rgb(62, 243, 227)';
	this.colors[87] = 'rgb(60, 245, 223)';
	this.colors[88] = 'rgb(58, 248, 219)';
	this.colors[89] = 'rgb(55, 251, 215)';
	this.colors[90] = 'rgb(53, 252, 211)';
	this.colors[91] = 'rgb(50, 253, 206)';
	this.colors[92] = 'rgb(47, 253, 201)';
	this.colors[93] = 'rgb(44, 253, 195)';
	this.colors[94] = 'rgb(42, 254, 190)';
	this.colors[95] = 'rgb(39, 254, 185)';
	this.colors[96] = 'rgb(36, 254, 180)';
	this.colors[97] = 'rgb(33, 254, 174)';
	this.colors[98] = 'rgb(31, 254, 168)';
	this.colors[99] = 'rgb(29, 253, 162)';
	this.colors[100] = 'rgb(26, 253, 156)';
	this.colors[101] = 'rgb(24, 252, 150)';
	this.colors[102] = 'rgb(22, 252, 144)';
	this.colors[103] = 'rgb(19, 252, 138)';
	this.colors[104] = 'rgb(18, 250, 131)';
	this.colors[105] = 'rgb(16, 248, 124)';
	this.colors[106] = 'rgb(15, 246, 118)';
	this.colors[107] = 'rgb(14, 243, 111)';
	this.colors[108] = 'rgb(12, 241, 104)';
	this.colors[109] = 'rgb(11, 239, 98)';
	this.colors[110] = 'rgb(10, 236, 91)';
	this.colors[111] = 'rgb(9, 234, 85)';
	this.colors[112] = 'rgb(8, 231, 78)';
	this.colors[113] = 'rgb(8, 229, 72)';
	this.colors[114] = 'rgb(7, 226, 66)';
	this.colors[115] = 'rgb(7, 224, 60)';
	this.colors[116] = 'rgb(6, 221, 54)';
	this.colors[117] = 'rgb(6, 219, 48)';
	this.colors[118] = 'rgb(7, 217, 43)';
	this.colors[119] = 'rgb(7, 216, 38)';
	this.colors[120] = 'rgb(8, 214, 33)';
	this.colors[121] = 'rgb(9, 213, 29)';
	this.colors[122] = 'rgb(10, 211, 24)';
	this.colors[123] = 'rgb(10, 210, 19)';
	this.colors[124] = 'rgb(11, 208, 15)';
	this.colors[125] = 'rgb(14, 209, 12)';
	this.colors[126] = 'rgb(16, 210, 10)';
	this.colors[127] = 'rgb(19, 211, 8)';
	this.colors[128] = 'rgb(21, 212, 6)';
	this.colors[129] = 'rgb(24, 213, 4)';
	this.colors[130] = 'rgb(27, 214, 2)';
	this.colors[131] = 'rgb(29, 215, 0)';
	this.colors[132] = 'rgb(34, 217, 0)';
	this.colors[133] = 'rgb(39, 220, 0)';
	this.colors[134] = 'rgb(44, 222, 0)';
	this.colors[135] = 'rgb(48, 224, 0)';
	this.colors[136] = 'rgb(53, 227, 0)';
	this.colors[137] = 'rgb(58, 229, 0)';
	this.colors[138] = 'rgb(63, 231, 0)';
	this.colors[139] = 'rgb(69, 234, 0)';
	this.colors[140] = 'rgb(74, 237, 0)';
	this.colors[141] = 'rgb(80, 239, 0)';
	this.colors[142] = 'rgb(86, 242, 0)';
	this.colors[143] = 'rgb(92, 245, 0)';
	this.colors[144] = 'rgb(98, 247, 0)';
	this.colors[145] = 'rgb(104, 250, 0)';
	this.colors[146] = 'rgb(111, 250, 0)';
	this.colors[147] = 'rgb(118, 251, 0)';
	this.colors[148] = 'rgb(125, 252, 0)';
	this.colors[149] = 'rgb(132, 253, 0)';
	this.colors[150] = 'rgb(139, 253, 0)';
	this.colors[151] = 'rgb(146, 254, 0)';
	this.colors[152] = 'rgb(153, 255, 0)';
	this.colors[153] = 'rgb(159, 255, 0)';
	this.colors[154] = 'rgb(165, 255, 0)';
	this.colors[155] = 'rgb(171, 255, 0)';
	this.colors[156] = 'rgb(177, 255, 0)';
	this.colors[157] = 'rgb(183, 255, 0)';
	this.colors[158] = 'rgb(189, 255, 0)';
	this.colors[159] = 'rgb(195, 255, 0)';
	this.colors[160] = 'rgb(200, 255, 0)';
	this.colors[161] = 'rgb(204, 255, 0)';
	this.colors[162] = 'rgb(209, 255, 0)';
	this.colors[163] = 'rgb(214, 255, 0)';
	this.colors[164] = 'rgb(219, 255, 0)';
	this.colors[165] = 'rgb(224, 255, 0)';
	this.colors[166] = 'rgb(227, 254, 0)';
	this.colors[167] = 'rgb(230, 252, 0)';
	this.colors[168] = 'rgb(233, 250, 0)';
	this.colors[169] = 'rgb(235, 249, 0)';
	this.colors[170] = 'rgb(238, 247, 0)';
	this.colors[171] = 'rgb(241, 245, 0)';
	this.colors[172] = 'rgb(243, 244, 0)';
	this.colors[173] = 'rgb(245, 240, 0)';
	this.colors[174] = 'rgb(247, 236, 0)';
	this.colors[175] = 'rgb(248, 232, 0)';
	this.colors[176] = 'rgb(250, 228, 0)';
	this.colors[177] = 'rgb(251, 224, 1)';
	this.colors[178] = 'rgb(253, 220, 1)';
	this.colors[179] = 'rgb(254, 216, 1)';
	this.colors[180] = 'rgb(254, 211, 1)';
	this.colors[181] = 'rgb(254, 206, 1)';
	this.colors[182] = 'rgb(254, 200, 1)';
	this.colors[183] = 'rgb(254, 195, 1)';
	this.colors[184] = 'rgb(254, 189, 1)';
	this.colors[185] = 'rgb(254, 184, 1)';
	this.colors[186] = 'rgb(254, 179, 2)';
	this.colors[187] = 'rgb(255, 173, 1)';
	this.colors[188] = 'rgb(255, 168, 1)';
	this.colors[189] = 'rgb(255, 162, 1)';
	this.colors[190] = 'rgb(255, 157, 0)';
	this.colors[191] = 'rgb(255, 151, 0)';
	this.colors[192] = 'rgb(255, 145, 0)';
	this.colors[193] = 'rgb(254, 140, 0)';
	this.colors[194] = 'rgb(254, 134, 0)';
	this.colors[195] = 'rgb(253, 129, 0)';
	this.colors[196] = 'rgb(252, 124, 0)';
	this.colors[197] = 'rgb(251, 118, 0)';
	this.colors[198] = 'rgb(250, 113, 0)';
	this.colors[199] = 'rgb(249, 108, 0)';
	this.colors[200] = 'rgb(248, 102, 0)';
	this.colors[201] = 'rgb(247, 98, 0)';
	this.colors[202] = 'rgb(245, 93, 0)';
	this.colors[203] = 'rgb(244, 88, 0)';
	this.colors[204] = 'rgb(243, 83, 0)';
	this.colors[205] = 'rgb(241, 79, 0)';
	this.colors[206] = 'rgb(240, 74, 0)';
	this.colors[207] = 'rgb(239, 70, 0)';
	this.colors[208] = 'rgb(237, 66, 0)';
	this.colors[209] = 'rgb(235, 63, 0)';
	this.colors[210] = 'rgb(234, 59, 0)';
	this.colors[211] = 'rgb(232, 56, 0)';
	this.colors[212] = 'rgb(230, 53, 0)';
	this.colors[213] = 'rgb(229, 49, 0)';
	this.colors[214] = 'rgb(227, 46, 0)';
	this.colors[215] = 'rgb(225, 45, 2)';
	this.colors[216] = 'rgb(224, 43, 4)';
	this.colors[217] = 'rgb(222, 42, 5)';
	this.colors[218] = 'rgb(221, 40, 7)';
	this.colors[219] = 'rgb(219, 38, 9)';
	this.colors[220] = 'rgb(217, 37, 11)';
	this.colors[221] = 'rgb(216, 36, 13)';
	this.colors[222] = 'rgb(215, 37, 16)';
	this.colors[223] = 'rgb(215, 38, 20)';
	this.colors[224] = 'rgb(214, 38, 23)';
	this.colors[225] = 'rgb(213, 39, 26)';
	this.colors[226] = 'rgb(212, 40, 29)';
	this.colors[227] = 'rgb(212, 40, 33)';
	this.colors[228] = 'rgb(212, 42, 36)';
	this.colors[229] = 'rgb(212, 46, 41)';
	this.colors[230] = 'rgb(213, 49, 45)';
	this.colors[231] = 'rgb(213, 52, 49)';
	this.colors[232] = 'rgb(214, 56, 53)';
	this.colors[233] = 'rgb(214, 59, 57)';
	this.colors[234] = 'rgb(215, 62, 62)';
	this.colors[235] = 'rgb(216, 67, 66)';
	this.colors[236] = 'rgb(218, 73, 71)';
	this.colors[237] = 'rgb(219, 79, 76)';
	this.colors[238] = 'rgb(221, 85, 81)';
	this.colors[239] = 'rgb(223, 90, 86)';
	this.colors[240] = 'rgb(224, 96, 91)';
	this.colors[241] = 'rgb(226, 102, 96)';
	this.colors[242] = 'rgb(228, 108, 100)';
	this.colors[243] = 'rgb(231, 115, 104)';
	this.colors[244] = 'rgb(233, 121, 107)';
	this.colors[245] = 'rgb(235, 128, 111)';
	this.colors[246] = 'rgb(238, 134, 115)';
	this.colors[247] = 'rgb(240, 141, 118)';
	this.colors[248] = 'rgb(243, 147, 122)';
	this.colors[249] = 'rgb(245, 153, 124)';
	this.colors[250] = 'rgb(246, 159, 125)';
	this.colors[251] = 'rgb(248, 165, 126)';
	this.colors[252] = 'rgb(250, 171, 128)';
	this.colors[253] = 'rgb(251, 177, 129)';
	this.colors[254] = 'rgb(253, 183, 130)';
	this.colors[255] = 'rgb(255, 189, 132)';
	if (this.projection) {
		this.startOffsetX = this.projection.offsetX;
		this.startOffsetY = this.projection.offsetY;
		this.startScale = this.projection.scale;
	}
};

MotionDisplay.prototype.setAlpha = function(alpha) {
	this.backgroundAlpha = 'rgba(' + this.rgb + ', ' + alpha + ')';
};

MotionDisplay.prototype.makeNewParticles = function(animator) {
	this.particles = [];
	for (var i = 0; i < this.numParticles; i++) {
		this.particles.push(this.makeParticle(animator));
	}
};

MotionDisplay.prototype.makeParticle = function(animator) {
	var dx = animator ? animator.dx : 0;
	var dy = animator ? animator.dy : 0;
	var scale = animator ? animator.scale : 1;
	var safecount = 0;
	for (;;) {
		var a = Math.random();
		var b = Math.random();
		var x = a * this.x0 + (1 - a) * this.x1;
		var y = b * this.y0 + (1 - b) * this.y1;
		var v = this.field.getValue(x, y);
		if (this.field.maxLength == 0) {
			return new Particle(x, y, 1 + 40 * Math.random());
		}
		var m = v.length() / this.field.maxLength;

		if ((v.x || v.y) && (++safecount > 10 || Math.random() > m * 0.9)) {
			var proj = this.projection.project(x, y);
			var sx = proj.x * scale + dx;
			var sy = proj.y * scale + dy;
			if (++safecount > 10 || !(sx < 0 || sy < 0 || sx > this.canvas.width || sy > this.canvas.height)) {
				return new Particle(x, y, 1 + 40 * Math.random());
			}
		}
	}
};

MotionDisplay.prototype.startMove = function(animator) {
	this.imageCanvas.getContext('2d').drawImage(this.canvas, 0, 0);
};

MotionDisplay.prototype.endMove = function(animator) {
	if (animator.scale < 1.1) {
		this.x0 = this.field.x0;
		this.x1 = this.field.x1;
		this.y0 = this.field.y0;
		this.y1 = this.field.y1;
	}
	else {
		var p = this.projection;
		var self = this;

		function invert(x, y) {
			x = (x - animator.dx) / animator.scale;
			y = (y - animator.dy) / animator.scale;
			return self.projection.invert(x, y);
		}

		var loc = invert(0, 0);
		var x0 = loc.x;
		var x1 = loc.x;
		var y0 = loc.y;
		var y1 = loc.y;

		function expand(x, y) {
			var v = invert(x, y);
			x0 = Math.min(v.x, x0);
			x1 = Math.max(v.x, x1);
			y0 = Math.min(v.y, y0);
			y1 = Math.max(v.y, y1);
		}

		var top = -0.2 * this.canvas.height;
		expand(top, this.canvas.height);
		expand(this.canvas.width, top);
		expand(this.canvas.width, this.canvas.height);
		this.x0 = Math.max(this.field.x0, x0);
		this.x1 = Math.min(this.field.x1, x1);
		this.y0 = Math.max(this.field.y0, y0);
		this.y1 = Math.min(this.field.y1, y1);
	}
	tick = 0;
	this.makeNewParticles(animator);
};

MotionDisplay.prototype.animate = function(animator) {
	this.moveThings(animator);
	this.draw(animator);
};

MotionDisplay.prototype.move = function(animator) {
	var w = this.canvas.width;
	var h = this.canvas.height;
	var g = this.canvas.getContext('2d');

	g.fillStyle = this.outsideColor;
	var dx = animator.dx;
	var dy = animator.dy;
	var scale = animator.scale;

	g.fillRect(0, 0, w, h);
	g.fillStyle = this.background;
	g.fillRect(dx, dy, w * scale, h * scale);
	var z = animator.relativeZoom();
	var dx = animator.dx - z * animator.dxStart;
	var dy = animator.dy - z * animator.dyStart;
	g.drawImage(this.imageCanvas, dx, dy, z * w, z * h);
};

MotionDisplay.prototype.moveThings = function(animator) {
	var speed = 0.01 * this.speedScale / animator.scale;
	for (var i = 0; i < this.particles.length; i++) {
		var p = this.particles[i];
		if (p.age > 0 && this.field.inBounds(p.x, p.y)) {
			var a = this.field.getValue(p.x, p.y);
			p.x += speed * a.x;
			p.y += speed * a.y;
			p.age--;
		}
		else {
			this.particles[i] = this.makeParticle(animator);
		}
	}
};

MotionDisplay.prototype.draw = function(animator) {
	var g = this.canvas.getContext('2d');
	var w = this.canvas.width;
	var h = this.canvas.height;
	if (this.first) {
		g.fillStyle = this.background;
		this.first = false;
	}
	else {
		g.fillStyle = this.backgroundAlpha;
	}
	var dx = animator.dx;
	var dy = animator.dy;
	var scale = animator.scale;

	g.fillRect(dx, dy, w * scale, h * scale);
	var proj = new Vector(0, 0);
	var val = new Vector(0, 0);
	g.lineWidth = 0.75;
	for (var i = 0; i < this.particles.length; i++) {
		var p = this.particles[i];
		if (!this.field.inBounds(p.x, p.y)) {
			p.age = -2;
			continue;
		}
		this.projection.project(p.x, p.y, proj);
		proj.x = proj.x * scale + dx;
		proj.y = proj.y * scale + dy;
		if (proj.x < 0 || proj.y < 0 || proj.x > w || proj.y > h) {
			p.age = -2;
		}
		if (p.oldX != -1) {
			var wind = this.field.getValue(p.x, p.y, val);
			var s = wind.length() / this.maxLength;
			var c = 90 + Math.round(70 * s);
			if (c > 255) {
				c = 255;
			}
			g.strokeStyle = this.colors[c];
			g.beginPath();
			g.moveTo(proj.x, proj.y);
			g.lineTo(p.oldX, p.oldY);
			g.stroke();
		}
		p.oldX = proj.x;
		p.oldY = proj.y;
	}
};

var MotionDetails = function(div, callout, field, projection, animator) {
	$(callout).fadeOut();
	var moveTime = +new Date();
	var calloutOK = false;
	var currentlyShowing = false;
	var calloutX = 0;
	var calloutY = 0;
	var calloutHTML = '';
	var lastX = 0;
	var lastY = 0;

	function format(x) {
		x = Math.round(x * 10) / 10;
		var a1 = ~~x;
		var a2 = (~~(x * 10)) % 10;
		return a1 + '.' + a2;
	}

	function minutes(x) {
		x = Math.round(x * 60) / 60;
		var degrees = ~~x;
		var m = ~~((x - degrees) * 60);
		return degrees + '&deg;&nbsp;' + (m == 0 ? '00' : m < 10 ? '0' + m : '' + m) + "'";
	}

	$(div).mouseleave(function() {
		moveTime = +new Date();
		calloutOK = false;
	});

	var pos = $(div).position();

	$(div).mousemove(function(e) {
		var x = e.pageX - this.offsetLeft - 60;
		var y = e.pageY - this.offsetTop - 10;
		if (x == lastX && y == lastY) {
			return;
		}
		lastX = x;
		lastY = y;
		moveTime = +new Date();
		var scale = animator.scale;
		var dx = animator.dx;
		var dy = animator.dy;
		var mx = (x - dx) / scale;
		var my = (y - dy) / scale;
		var location = projection.invert(mx, my);
		var lat = location.y;
		var lon = location.x;
		var speed = 0;
		if (field.inBounds(lon, lat)) {
			speed = field.getValue(lon, lat).length() * 2.23693629;
		}
		calloutOK = !!speed;
		calloutHTML = 	'<div style="padding-bottom:5px"><b>' +
						format(speed) + ' mph</b> wind speed<br></div>' +
						minutes(lat) + ' N, ' +
						minutes(-lon) + ' W<br>' +
						'click to zoom';
		calloutY = (pos.top + y) + 'px';
		calloutX = (pos.left + x + 20) + 'px';
	});

	setInterval(function() {
		var timeSinceMove = +new Date() - moveTime;
		if (timeSinceMove > 200 && calloutOK) {
			if (!currentlyShowing) {
				callout.innerHTML = calloutHTML;
				callout.style.left = calloutX;
				callout.style.top = calloutY;
				callout.style.visibility = 'visible';
				$(callout).fadeTo(400, 1);
				currentlyShowing = true;
			}
		}
		else if (currentlyShowing) {
			$(callout).fadeOut('fast');
			currentlyShowing = false;
		}
	}, 50);
};

var CityDisplay = function(cities, canvas, projection) {
	this.cities = cities;
	this.canvas = canvas;
	this.projection = projection;
	this.maxInView = 10;
	this.pad = 3;
	cities.sort(function(a, b) {
		return b.pop - a.pop;
	});
	for (var i = 0; i < this.cities.length; i++) {
		this.cities[i].alpha = 0;
	}
};

CityDisplay.prototype.endMove = function(animator) {
	for (var i = 0; i < this.cities.length; i++) {
		this.cities[i].alpha = 0;
	}
	this.move(animator);
};

CityDisplay.prototype.markCities = function(scale, dx, dy, alpha) {
	var spaceTaken = [];
	function collide(r1, r2) {
		return	!(r1.x + r1.w < r2.x || r1.x > r2.x + r2.w ||
				  r1.y + r1.h < r2.y || r1.y > r2.y + r2.h);
	}

	function isFree(r) {
		for (var i = 0; i < spaceTaken.length; i++) {
			if (collide(r, spaceTaken[i])) {
				return false;
			}
		}
		return true;
	}

	var g = this.canvas.getContext('2d');
	var w = this.canvas.width;
	var h = this.canvas.height;
	var numInView = 0;
	var pad = this.pad;
	for (var i = 0; i < this.cities.length; i++) {
		var city = this.cities[i];
		var r = 0.075 * Math.pow(city.pop, 0.3);
		var v = this.projection.project(city.lon, city.lat);
		var x = v.x * scale + dx;
		var y = v.y * scale + dy;
		
		if (x < 0 || x > w || y < 0 || y > h) {
			continue;
		}

		var tx = x;
		var ty = y + 15;
		var textSize = g.measureText(city.city);
		
		var dotArea = {
			'x': x - r - pad,
			'y': y - r - pad,
			'w': 2 * (r + pad),
			'h': 2 * (r + pad)
		};
		
		var textArea = {
			'x': tx - textSize.width / 2 - pad,
			'y': ty - 15 - pad,
			'w': textSize.width + 2 * pad,
			'h': 15 + 2 * pad
		};

		if (!isFree(dotArea) || !isFree(textArea)) {
			continue;
		}

		spaceTaken.push(textArea);
		spaceTaken.push(dotArea);
		city.alpha += alpha;
		if (++numInView > this.maxInView) {
			break;
		}
	}
};

CityDisplay.prototype.move = function(animator) {
	for (var i = 0; i < this.cities.length; i++) {
		this.cities[i].alpha = 0;
	}
	var dx = 0;
	var dy = 0;
	var scale = 1;
	
	if (animator) {
		dx = animator.dx;
		dy = animator.dy;
		scale = animator.scale;
		if (animator.state == 'zoom') {
			var u = animator.zoomProgress;
			this.markCities(animator.scaleStart, animator.dxStart, animator.dyStart, 1 - u);
			this.markCities(animator.scaleTarget, animator.dxTarget, animator.dyTarget, u);
		}
		else {
			this.markCities(scale, dx, dy, 1);
		}
	}
	else {
		this.markCities(1, 0, 0, 1);
	}

	var g = this.canvas.getContext('2d');
	var w = this.canvas.width;
	var h = this.canvas.height;
	g.clearRect(0, 0, w, h);
	var pad = this.pad;
	for (var i = 0; i < this.cities.length; i++) {
		var city = this.cities[i];
		var alpha = Math.min(1, city.alpha);
		if (!alpha) {
			continue;
		}

		function check(val, name) {
			if (!val) {
				window.console.log(name + ' = ' + val);
			}
		}

		var r = 0.075 * Math.pow(city.pop, 0.3);
		var v = this.projection.project(city.lon, city.lat);
		var x = v.x * scale + dx;
		var y = v.y * scale + dy;

		if (x < 0 || x > w || y < 0 || y > h) {
			continue;
		}

		var tx = x;
		var ty = y + 15;

		g.beginPath();
		g.arc(x, y, r, 0, Math.PI * 2, true);
		g.closePath();
		g.fillStyle = 'rgba(255, 255, 255, ' + alpha + ')';
		g.fill();
		g.strokeStyle = 'rgba(0, 0, 0, ' + alpha + ')';
		g.stroke();
		g.fillStyle = 'rgba(255, 255, 255, ' + 0.25 * alpha + ')';
		g.textAlign = 'center';
		g.font = '12px Verdana';

		for (var a = -2; a <= 2; a++) {
			for (var b = -2; b <= 2; b++) {
				g.fillText(city.city, tx + a, ty + b);
			}
		}
		g.fillStyle = 'rgba(0, 0, 0, ' + alpha + ')';
		g.fillText(city.city, tx, ty);
	}
};

var cities = [
	{
		city: 'Abbyville',
		state: 'Kansas',
		lat: 37.9702735,
		lon: -98.204054584378525,
		pop: 87
	},
	{
		city: 'Abilene',
		state: 'Kansas',
		lat: 38.9261675,
		lon: -97.213930257507954,
		pop: 6844
	},
	{
		city: 'Ada',
		state: 'Kansas',
		lat: 39.160951,
		lon: -97.886908023854161,
		pop: 100
	},
	{
		city: 'Admire',
		state: 'Kansas',
		lat: 38.6425345,
		lon: -96.102068917388848,
		pop: 156
	},
	{
		city: 'Agenda',
		state: 'Kansas',
		lat: 39.7059755,
		lon: -97.431962082335303,
		pop: 68
	},
	{
		city: 'Agra',
		state: 'Kansas',
		lat: 39.761103,
		lon: -99.119663304904833,
		pop: 267
	},
	{
		city: 'Albert',
		state: 'Kansas',
		lat: 38.4539585,
		lon: -99.011025874494749,
		pop: 175
	},
	{
		city: 'Alden',
		state: 'Kansas',
		lat: 38.2411425,
		lon: -98.312014655155522,
		pop: 148
	},
	{
		city: 'Alexander',
		state: 'Kansas',
		lat: 38.4692475,
		lon: -99.552466452189933,
		pop: 65
	},
	{
		city: 'Allen',
		state: 'Kansas',
		lat: 38.6546555,
		lon: -96.169667335167318,
		pop: 177
	},
	{
		city: 'Alma',
		state: 'Kansas',
		lat: 39.01477,
		lon: -96.288831446395562,
		pop: 832
	},
	{
		city: 'Almena',
		state: 'Kansas',
		lat: 39.8911025,
		lon: -99.709591585923704,
		pop: 408
	},
	{
		city: 'Alta Vista',
		state: 'Kansas',
		lat: 38.862342,
		lon: -96.488148363119066,
		pop: 444
	},
	{
		city: 'Altamont',
		state: 'Kansas',
		lat: 37.1933805,
		lon: -95.29429741226464,
		pop: 1080
	},
	{
		city: 'Alton',
		state: 'Kansas',
		lat: 39.4686195,
		lon: -98.948102196874572,
		pop: 103
	},
	{
		city: 'Altoona',
		state: 'Kansas',
		lat: 37.526938,
		lon: -95.661158881722685,
		pop: 414
	},
	{
		city: 'Americus',
		state: 'Kansas',
		lat: 38.5049105,
		lon: -96.26178075,
		pop: 894
	},
	{
		city: 'Andale',
		state: 'Kansas',
		lat: 37.7917555,
		lon: -97.625366446174851,
		pop: 928
	},
	{
		city: 'Andover',
		state: 'Kansas',
		lat: 37.683234,
		lon: -97.139527756186823,
		pop: 11791
	},
	{
		city: 'Anthony',
		state: 'Kansas',
		lat: 37.1552075,
		lon: -98.03242503765324,
		pop: 2269
	},
	{
		city: 'Arcadia',
		state: 'Kansas',
		lat: 37.642306,
		lon: -94.624044922446402,
		pop: 310
	},
	{
		city: 'Argonia',
		state: 'Kansas',
		lat: 37.2669165,
		lon: -97.765126420634914,
		pop: 501
	},
	{
		city: 'Arkansas City',
		state: 'Kansas',
		lat: 37.076312,
		lon: -97.037152011440128,
		pop: 12415
	},
	{
		city: 'Arlington',
		state: 'Kansas',
		lat: 37.895446,
		lon: -98.177444806947477,
		pop: 473
	},
	{
		city: 'Arma',
		state: 'Kansas',
		lat: 37.5408105,
		lon: -94.700354414846259,
		pop: 1481
	},
	{
		city: 'Asherville',
		state: 'Kansas',
		lat: 39.407444,
		lon: -97.972001376425411,
		pop: 28
	},
	{
		city: 'Ashland',
		state: 'Kansas',
		lat: 37.1861895,
		lon: -99.76731841933497,
		pop: 867
	},
	{
		city: 'Assaria',
		state: 'Kansas',
		lat: 38.680086,
		lon: -97.604458181818188,
		pop: 413
	},
	{
		city: 'Atchison',
		state: 'Kansas',
		lat: 39.561281,
		lon: -95.133064817761991,
		pop: 11021
	},
	{
		city: 'Athol',
		state: 'Kansas',
		lat: 39.7654495,
		lon: -98.919993075574808,
		pop: 44
	},
	{
		city: 'Atlanta',
		state: 'Kansas',
		lat: 37.4349775,
		lon: -96.767221104092002,
		pop: 195
	},
	{
		city: 'Attica',
		state: 'Kansas',
		lat: 37.243499,
		lon: -98.226039675,
		pop: 626
	},
	{
		city: 'Atwood',
		state: 'Kansas',
		lat: 39.8118495,
		lon: -101.040563773327136,
		pop: 1194
	},
	{
		city: 'Auburn',
		state: 'Kansas',
		lat: 38.9073175,
		lon: -95.816945229618526,
		pop: 1227
	},
	{
		city: 'Augusta',
		state: 'Kansas',
		lat: 37.696181,
		lon: -96.973271628889165,
		pop: 9274
	},
	{
		city: 'Aurora',
		state: 'Kansas',
		lat: 39.452047,
		lon: -97.53033871,
		pop: 60
	},
	{
		city: 'Axtell',
		state: 'Kansas',
		lat: 39.8703755,
		lon: -96.254737229775657,
		pop: 406
	},
	{
		city: 'Baileyville',
		state: 'Kansas',
		lat: 39.8419245,
		lon: -96.182872765373787,
		pop: 181
	},
	{
		city: 'Baldwin City',
		state: 'Kansas',
		lat: 38.781275,
		lon: -95.19053472158528,
		pop: 4515
	},
	{
		city: 'Barnard',
		state: 'Kansas',
		lat: 39.189073,
		lon: -98.043526993945648,
		pop: 70
	},
	{
		city: 'Barnes',
		state: 'Kansas',
		lat: 39.712265,
		lon: -96.873095386016658,
		pop: 159
	},
	{
		city: 'Bartlett',
		state: 'Kansas',
		lat: 37.055077,
		lon: -95.211108413071244,
		pop: 80
	},
	{
		city: 'Basehor',
		state: 'Kansas',
		lat: 39.1298555,
		lon: -94.928894082223138,
		pop: 4613
	},
	{
		city: 'Bassett',
		state: 'Kansas',
		lat: 37.9058715,
		lon: -95.407511004581153,
		pop: 14
	},
	{
		city: 'Baxter Springs',
		state: 'Kansas',
		lat: 37.020304,
		lon: -94.735550949387729,
		pop: 4238
	},
	{
		city: 'Bazine',
		state: 'Kansas',
		lat: 38.4455975,
		lon: -99.693083577313672,
		pop: 334
	},
	{
		city: 'Beattie',
		state: 'Kansas',
		lat: 39.862056,
		lon: -96.417856827285419,
		pop: 200
	},
	{
		city: 'Bel Aire',
		state: 'Kansas',
		lat: 37.773995,
		lon: -97.253670592107397,
		pop: 6769
	},
	{
		city: 'Belle Plaine',
		state: 'Kansas',
		lat: 37.3955635,
		lon: -97.279637228458483,
		pop: 1681
	},
	{
		city: 'Belleville',
		state: 'Kansas',
		lat: 39.825381,
		lon: -97.626187497146489,
		pop: 1991
	},
	{
		city: 'Beloit',
		state: 'Kansas',
		lat: 39.4649995,
		lon: -98.106597129537604,
		pop: 3835
	},
	{
		city: 'Belpre',
		state: 'Kansas',
		lat: 37.951138,
		lon: -99.09924072317267,
		pop: 84
	},
	{
		city: 'Belvue',
		state: 'Kansas',
		lat: 39.216939,
		lon: -96.177544311127704,
		pop: 205
	},
	{
		city: 'Bendena',
		state: 'Kansas',
		lat: 39.7458895,
		lon: -95.181030094578475,
		pop: 117
	},
	{
		city: 'Benedict',
		state: 'Kansas',
		lat: 37.627533,
		lon: -95.744552110884342,
		pop: 73
	},
	{
		city: 'Bennington',
		state: 'Kansas',
		lat: 39.0337545,
		lon: -97.593704581030991,
		pop: 672
	},
	{
		city: 'Bentley',
		state: 'Kansas',
		lat: 37.886802,
		lon: -97.516614778163415,
		pop: 530
	},
	{
		city: 'Benton',
		state: 'Kansas',
		lat: 37.7853285,
		lon: -97.107911914119683,
		pop: 880
	},
	{
		city: 'Bern',
		state: 'Kansas',
		lat: 39.9609555,
		lon: -95.971677291959736,
		pop: 166
	},
	{
		city: 'Beverly',
		state: 'Kansas',
		lat: 39.0135515,
		lon: -97.974741000969331,
		pop: 162
	},
	{
		city: 'Bird City',
		state: 'Kansas',
		lat: 39.746227,
		lon: -101.53260804436934,
		pop: 447
	},
	{
		city: 'Bison',
		state: 'Kansas',
		lat: 38.5194025,
		lon: -99.19769096133416,
		pop: 255
	},
	{
		city: 'Blue Mound',
		state: 'Kansas',
		lat: 38.0882955,
		lon: -95.009384726682896,
		pop: 275
	},
	{
		city: 'Blue Rapids',
		state: 'Kansas',
		lat: 39.680184,
		lon: -96.659949167832167,
		pop: 1019
	},
	{
		city: 'Bluff City',
		state: 'Kansas',
		lat: 37.076952,
		lon: -97.874947290218628,
		pop: 65
	},
	{
		city: 'Bogue',
		state: 'Kansas',
		lat: 39.360606,
		lon: -99.68748987068966,
		pop: 143
	},
	{
		city: 'Bonner Springs',
		state: 'Kansas',
		lat: 39.077812,
		lon: -94.871842244703174,
		pop: 7314
	},
	{
		city: 'Brewster',
		state: 'Kansas',
		lat: 39.362667,
		lon: -101.377121030651935,
		pop: 305
	},
	{
		city: 'Bronson',
		state: 'Kansas',
		lat: 37.897246,
		lon: -95.072491413501595,
		pop: 323
	},
	{
		city: 'Brookville',
		state: 'Kansas',
		lat: 38.774264,
		lon: -97.86409531368821,
		pop: 262
	},
	{
		city: 'Brownell',
		state: 'Kansas',
		lat: 38.6400545,
		lon: -99.745018743998884,
		pop: 29
	},
	{
		city: 'Bucklin',
		state: 'Kansas',
		lat: 37.549297,
		lon: -99.636009628378375,
		pop: 794
	},
	{
		city: 'Bucyrus',
		state: 'Kansas',
		lat: 38.7237035,
		lon: -94.713993663451816,
		pop: 193
	},
	{
		city: 'Buffalo',
		state: 'Kansas',
		lat: 37.7086835,
		lon: -95.697327286647308,
		pop: 232
	},
	{
		city: 'Buhler',
		state: 'Kansas',
		lat: 38.139865,
		lon: -97.769278712802418,
		pop: 1327
	},
	{
		city: 'Bunker Hill',
		state: 'Kansas',
		lat: 38.8760455,
		lon: -98.699568028447928,
		pop: 95
	},
	{
		city: 'Burden',
		state: 'Kansas',
		lat: 37.314478,
		lon: -96.756089703821658,
		pop: 535
	},
	{
		city: 'Burdett',
		state: 'Kansas',
		lat: 38.1935015,
		lon: -99.526987238556174,
		pop: 247
	},
	{
		city: 'Burlingame',
		state: 'Kansas',
		lat: 38.751113,
		lon: -95.834204863193165,
		pop: 934
	},
	{
		city: 'Burlington',
		state: 'Kansas',
		lat: 38.1941865,
		lon: -95.74594410248433,
		pop: 2674
	},
	{
		city: 'Burns',
		state: 'Kansas',
		lat: 38.088026,
		lon: -96.887899113229011,
		pop: 228
	},
	{
		city: 'Burr Oak',
		state: 'Kansas',
		lat: 39.869628,
		lon: -98.305999574497434,
		pop: 174
	},
	{
		city: 'Burrton',
		state: 'Kansas',
		lat: 38.024971,
		lon: -97.671868414634147,
		pop: 901
	},
	{
		city: 'Bushong',
		state: 'Kansas',
		lat: 38.6427085,
		lon: -96.256675387126862,
		pop: 34
	},
	{
		city: 'Bushton',
		state: 'Kansas',
		lat: 38.513002,
		lon: -98.39603195942064,
		pop: 279
	},
	{
		city: 'Byers',
		state: 'Kansas',
		lat: 37.787411,
		lon: -98.866975074934459,
		pop: 35
	},
	{
		city: 'Caldwell',
		state: 'Kansas',
		lat: 37.03432,
		lon: -97.609735415098655,
		pop: 1068
	},
	{
		city: 'Cambridge',
		state: 'Kansas',
		lat: 37.316988,
		lon: -96.666332246636756,
		pop: 82
	},
	{
		city: 'Caney',
		state: 'Kansas',
		lat: 37.0134935,
		lon: -95.932758882940462,
		pop: 2203
	},
	{
		city: 'Canton',
		state: 'Kansas',
		lat: 38.3848815,
		lon: -97.430281184234161,
		pop: 748
	},
	{
		city: 'Carbondale',
		state: 'Kansas',
		lat: 38.818975,
		lon: -95.694462652913018,
		pop: 1437
	},
	{
		city: 'Carlton',
		state: 'Kansas',
		lat: 38.686608,
		lon: -97.29319947687415,
		pop: 42
	},
	{
		city: 'Cassoday',
		state: 'Kansas',
		lat: 38.039944,
		lon: -96.637502682640147,
		pop: 129
	},
	{
		city: 'Catharine',
		state: 'Kansas',
		lat: 38.9292435,
		lon: -99.21488752231997,
		pop: 104
	},
	{
		city: 'Cawker City',
		state: 'Kansas',
		lat: 39.50988,
		lon: -98.433082680875572,
		pop: 469
	},
	{
		city: 'Cedar',
		state: 'Kansas',
		lat: 39.6577985,
		lon: -98.939977345469259,
		pop: 28
	},
	{
		city: 'Cedar Point',
		state: 'Kansas',
		lat: 38.2596755,
		lon: -96.821524105227383,
		pop: 579
	},
	{
		city: 'Cedar Vale',
		state: 'Kansas',
		lat: 37.1067,
		lon: -96.501680539201175,
		pop: 14
	},
	{
		city: 'Centralia',
		state: 'Kansas',
		lat: 39.7249375,
		lon: -96.129479511090523,
		pop: 512
	},
	{
		city: 'Chanute',
		state: 'Kansas',
		lat: 37.666734,
		lon: -95.458570136127307,
		pop: 9119
	},
	{
		city: 'Chapman',
		state: 'Kansas',
		lat: 38.9737985,
		lon: -97.024919969972089,
		pop: 1393
	},
	{
		city: 'Chase',
		state: 'Kansas',
		lat: 38.3550115,
		lon: -98.349299943217801,
		pop: 477
	},
	{
		city: 'Chautauqua',
		state: 'Kansas',
		lat: 37.0250905,
		lon: -96.17764632275177,
		pop: 111
	},
	{
		city: 'Cheney',
		state: 'Kansas',
		lat: 37.632886,
		lon: -97.780895900674864,
		pop: 2094
	},
	{
		city: 'Cherokee',
		state: 'Kansas',
		lat: 37.3463385,
		lon: -94.818533880202779,
		pop: 714
	},
	{
		city: 'Cherryvale',
		state: 'Kansas',
		lat: 37.268012,
		lon: -95.552109321818961,
		pop: 2367
	},
	{
		city: 'Chetopa',
		state: 'Kansas',
		lat: 37.0390915,
		lon: -95.095282662267095,
		pop: 1125
	},
	{
		city: 'Chicopee',
		state: 'Kansas',
		lat: 37.3826245,
		lon: -94.741753019477542,
		pop: 408
	},
	{
		city: 'Cimarron',
		state: 'Kansas',
		lat: 37.809816,
		lon: -100.345863201437993,
		pop: 2184
	},
	{
		city: 'Circleville',
		state: 'Kansas',
		lat: 39.5095935,
		lon: -95.855495510871407,
		pop: 170
	},
	{
		city: 'Claflin',
		state: 'Kansas',
		lat: 38.5250265,
		lon: -98.536596517270496,
		pop: 645
	},
	{
		city: 'Clay Center',
		state: 'Kansas',
		lat: 39.383061,
		lon: -97.119926257491031,
		pop: 4334
	},
	{
		city: 'Clayton',
		state: 'Kansas',
		lat: 39.738292,
		lon: -100.176730508740519,
		pop: 59
	},
	{
		city: 'Clearwater',
		state: 'Kansas',
		lat: 37.5082305,
		lon: -97.499083544642843,
		pop: 2481
	},
	{
		city: 'Clifton',
		state: 'Kansas',
		lat: 39.568353,
		lon: -97.279845376759866,
		pop: 554
	},
	{
		city: 'Climax',
		state: 'Kansas',
		lat: 37.7183705,
		lon: -96.223356999989619,
		pop: 72
	},
	{
		city: 'Clyde',
		state: 'Kansas',
		lat: 39.591367,
		lon: -97.401620010720123,
		pop: 716
	},
	{
		city: 'Coats',
		state: 'Kansas',
		lat: 37.511071,
		lon: -98.82497275,
		pop: 83
	},
	{
		city: 'Coffeyville',
		state: 'Kansas',
		lat: 37.040779,
		lon: -95.610154187120187,
		pop: 10295
	},
	{
		city: 'Colby',
		state: 'Kansas',
		lat: 39.380655,
		lon: -101.046514076142074,
		pop: 5387
	},
	{
		city: 'Coldwater',
		state: 'Kansas',
		lat: 37.257484,
		lon: -99.346429422071367,
		pop: 828
	},
	{
		city: 'Collyer',
		state: 'Kansas',
		lat: 39.036095,
		lon: -100.117989092155625,
		pop: 109
	},
	{
		city: 'Colony',
		state: 'Kansas',
		lat: 38.070656,
		lon: -95.363931008016237,
		pop: 408
	},
	{
		city: 'Columbus',
		state: 'Kansas',
		lat: 37.175406,
		lon: -94.838377772953024,
		pop: 3312
	},
	{
		city: 'Colwich',
		state: 'Kansas',
		lat: 37.7828295,
		lon: -97.536062423594359,
		pop: 1327
	},
	{
		city: 'Concordia',
		state: 'Kansas',
		lat: 39.561017,
		lon: -97.658864865720119,
		pop: 5395
	},
	{
		city: 'Conway Springs',
		state: 'Kansas',
		lat: 37.3893335,
		lon: -97.642952590336265,
		pop: 1272
	},
	{
		city: 'Coolidge',
		state: 'Kansas',
		lat: 38.0411105,
		lon: -102.008519257978719,
		pop: 95
	},
	{
		city: 'Copeland',
		state: 'Kansas',
		lat: 37.53934,
		lon: -100.629376628075931,
		pop: 310
	},
	{
		city: 'Corning',
		state: 'Kansas',
		lat: 39.6567355,
		lon: -96.028987004897488,
		pop: 157
	},
	{
		city: 'Cottonwood Falls',
		state: 'Kansas',
		lat: 38.3697245,
		lon: -96.541393059585488,
		pop: 903
	},
	{
		city: 'Council Grove',
		state: 'Kansas',
		lat: 38.6647325,
		lon: -96.489455913712163,
		pop: 2182
	},
	{
		city: 'Courtland',
		state: 'Kansas',
		lat: 39.7818655,
		lon: -97.895468277553405,
		pop: 285
	},
	{
		city: 'Coyville',
		state: 'Kansas',
		lat: 37.687792,
		lon: -95.893876037296735,
		pop: 46
	},
	{
		city: 'Cuba',
		state: 'Kansas',
		lat: 39.801234,
		lon: -97.457914113103783,
		pop: 156
	},
	{
		city: 'Cullison',
		state: 'Kansas',
		lat: 37.6293675,
		lon: -98.906208488483131,
		pop: 101
	},
	{
		city: 'Culver',
		state: 'Kansas',
		lat: 38.9682085,
		lon: -97.759012727471827,
		pop: 121
	},
	{
		city: 'Cunningham',
		state: 'Kansas',
		lat: 37.644403,
		lon: -98.43325141162029,
		pop: 454
	},
	{
		city: 'Damar',
		state: 'Kansas',
		lat: 39.318517,
		lon: -99.584566758093899,
		pop: 132
	},
	{
		city: 'Danville',
		state: 'Kansas',
		lat: 37.2866215,
		lon: -97.891849534812351,
		pop: 38
	},
	{
		city: 'De Soto',
		state: 'Kansas',
		lat: 38.9716405,
		lon: -94.955087356731525,
		pop: 5720
	},
	{
		city: 'Dearing',
		state: 'Kansas',
		lat: 37.0523565,
		lon: -95.69327060611144,
		pop: 431
	},
	{
		city: 'Deerfield',
		state: 'Kansas',
		lat: 37.982166,
		lon: -101.132145105925517,
		pop: 700
	},
	{
		city: 'Delia',
		state: 'Kansas',
		lat: 39.238922,
		lon: -95.965700952270907,
		pop: 169
	},
	{
		city: 'Delphos',
		state: 'Kansas',
		lat: 39.275236,
		lon: -97.766503489247299,
		pop: 359
	},
	{
		city: 'Denison',
		state: 'Kansas',
		lat: 39.3932175,
		lon: -95.62854603642981,
		pop: 187
	},
	{
		city: 'Denton',
		state: 'Kansas',
		lat: 39.7326565,
		lon: -95.271891090072188,
		pop: 148
	},
	{
		city: 'Derby',
		state: 'Kansas',
		lat: 37.555252,
		lon: -97.249057788939552,
		pop: 22158
	},
	{
		city: 'Detroit',
		state: 'Kansas',
		lat: 38.935239,
		lon: -97.123640304675263,
		pop: 114
	},
	{
		city: 'Dexter',
		state: 'Kansas',
		lat: 37.180019,
		lon: -96.715804264635096,
		pop: 278
	},
	{
		city: 'Dighton',
		state: 'Kansas',
		lat: 38.482058,
		lon: -100.466250783451073,
		pop: 1038
	},
	{
		city: 'Dodge City',
		state: 'Kansas',
		lat: 37.7659585,
		lon: -100.023909102661591,
		pop: 27340
	},
	{
		city: 'Dorrance',
		state: 'Kansas',
		lat: 38.848228,
		lon: -98.590515564438746,
		pop: 185
	},
	{
		city: 'Douglass',
		state: 'Kansas',
		lat: 37.5183965,
		lon: -97.00940121808344,
		pop: 1700
	},
	{
		city: 'Downs',
		state: 'Kansas',
		lat: 39.5020365,
		lon: -98.544708980769158,
		pop: 900
	},
	{
		city: 'Dresden',
		state: 'Kansas',
		lat: 39.62167,
		lon: -100.418924636960142,
		pop: 41
	},
	{
		city: 'Dunlap',
		state: 'Kansas',
		lat: 38.5761345,
		lon: -96.36758254123481,
		pop: 30
	},
	{
		city: 'Durham',
		state: 'Kansas',
		lat: 38.484965,
		lon: -97.227440754495177,
		pop: 112
	},
	{
		city: 'Dwight',
		state: 'Kansas',
		lat: 38.8447615,
		lon: -96.592287434664826,
		pop: 272
	},
	{
		city: 'Earlton',
		state: 'Kansas',
		lat: 37.5874245,
		lon: -95.469826431818177,
		pop: 55
	},
	{
		city: 'Eastborough',
		state: 'Kansas',
		lat: 37.6843775,
		lon: -97.257920611526743,
		pop: 773
	},
	{
		city: 'Easton',
		state: 'Kansas',
		lat: 39.344341,
		lon: -95.11604124316159,
		pop: 253
	},
	{
		city: 'Edgerton',
		state: 'Kansas',
		lat: 38.751348,
		lon: -95.001162696048127,
		pop: 1671
	},
	{
		city: 'Edmond',
		state: 'Kansas',
		lat: 39.626783,
		lon: -99.81939298040723,
		pop: 49
	},
	{
		city: 'Edna',
		state: 'Kansas',
		lat: 37.056498,
		lon: -95.359049999656165,
		pop: 442
	},
	{
		city: 'Edwardsville',
		state: 'Kansas',
		lat: 39.07477,
		lon: -94.816270099550422,
		pop: 4340
	},
	{
		city: 'Effingham',
		state: 'Kansas',
		lat: 39.5220875,
		lon: -95.398044342385518,
		pop: 546
	},
	{
		city: 'El Dorado',
		state: 'Kansas',
		lat: 37.8217565,
		lon: -96.870179741854628,
		pop: 13021
	},
	{
		city: 'Elbing',
		state: 'Kansas',
		lat: 38.0544915,
		lon: -97.127064849013919,
		pop: 229
	},
	{
		city: 'Elgin',
		state: 'Kansas',
		lat: 37.001812,
		lon: -96.280977574254749,
		pop: 89
	},
	{
		city: 'Elk City',
		state: 'Kansas',
		lat: 37.290263,
		lon: -95.910159141696028,
		pop: 325
	},
	{
		city: 'Elk Falls',
		state: 'Kansas',
		lat: 37.3741505,
		lon: -96.192877659158881,
		pop: 107
	},
	{
		city: 'Elkhart',
		state: 'Kansas',
		lat: 37.005104,
		lon: -101.894181881544625,
		pop: 2205
	},
	{
		city: 'Ellinwood',
		state: 'Kansas',
		lat: 38.3567015,
		lon: -98.581401758390669,
		pop: 2131
	},
	{
		city: 'Ellis',
		state: 'Kansas',
		lat: 38.935263,
		lon: -99.559640904422253,
		pop: 2062
	},
	{
		city: 'Ellsworth',
		state: 'Kansas',
		lat: 38.736654,
		lon: -98.227382005214224,
		pop: 3120
	},
	{
		city: 'Elmdale',
		state: 'Kansas',
		lat: 38.373811,
		lon: -96.645055614319133,
		pop: 55
	},
	{
		city: 'Elsmore',
		state: 'Kansas',
		lat: 37.7944035,
		lon: -95.149850717863075,
		pop: 77
	},
	{
		city: 'Elwood',
		state: 'Kansas',
		lat: 39.746944,
		lon: -94.887281909159185,
		pop: 1224
	},
	{
		city: 'Emmett',
		state: 'Kansas',
		lat: 39.30734,
		lon: -96.055655061053699,
		pop: 191
	},
	{
		city: 'Emporia',
		state: 'Kansas',
		lat: 38.411016,
		lon: -96.196099359130528,
		pop: 24916
	},
	{
		city: 'Englewood',
		state: 'Kansas',
		lat: 37.0395345,
		lon: -99.986090579514553,
		pop: 77
	},
	{
		city: 'Ensign',
		state: 'Kansas',
		lat: 37.653906,
		lon: -100.232912521137166,
		pop: 187
	},
	{
		city: 'Enterprise',
		state: 'Kansas',
		lat: 38.9021505,
		lon: -97.116446301549288,
		pop: 855
	},
	{
		city: 'Erie',
		state: 'Kansas',
		lat: 37.5739195,
		lon: -95.244459368284524,
		pop: 1150
	},
	{
		city: 'Esbon',
		state: 'Kansas',
		lat: 39.8224095,
		lon: -98.434356720806633,
		pop: 99
	},
	{
		city: 'Eskridge',
		state: 'Kansas',
		lat: 38.860306,
		lon: -96.105456049586024,
		pop: 534
	},
	{
		city: 'Eudora',
		state: 'Kansas',
		lat: 38.9376615,
		lon: -95.094819587409205,
		pop: 6136
	},
	{
		city: 'Eureka',
		state: 'Kansas',
		lat: 37.826206,
		lon: -96.288381712341192,
		pop: 2633
	},
	{
		city: 'Everest',
		state: 'Kansas',
		lat: 39.6764375,
		lon: -95.424639934820647,
		pop: 284
	},
	{
		city: 'Fairview',
		state: 'Kansas',
		lat: 39.840032,
		lon: -95.727493897435892,
		pop: 260
	},
	{
		city: 'Fairway',
		state: 'Kansas',
		lat: 39.0248775,
		lon: -94.630122627705191,
		pop: 3882
	},
	{
		city: 'Fall River',
		state: 'Kansas',
		lat: 37.6082245,
		lon: -96.028753863890529,
		pop: 162
	},
	{
		city: 'Falun',
		state: 'Kansas',
		lat: 38.6747765,
		lon: -97.751368227575398,
		pop: 87
	},
	{
		city: 'Florence',
		state: 'Kansas',
		lat: 38.244676,
		lon: -96.929168735529458,
		pop: 465
	},
	{
		city: 'Fontana',
		state: 'Kansas',
		lat: 38.422559,
		lon: -94.853179309148913,
		pop: 224
	},
	{
		city: 'Ford',
		state: 'Kansas',
		lat: 37.636962,
		lon: -99.75376004207331,
		pop: 216
	},
	{
		city: 'Formoso',
		state: 'Kansas',
		lat: 39.779145,
		lon: -97.993837819498069,
		pop: 93
	},
	{
		city: 'Fort Dodge',
		state: 'Kansas',
		lat: 37.7312425,
		lon: -99.937602770518154,
		pop: 165
	},
	{
		city: 'Fort Riley',
		state: 'Kansas',
		lat: 39.108705,
		lon: -96.814866762279593,
		pop: 7761
	},
	{
		city: 'Fort Scott',
		state: 'Kansas',
		lat: 37.827899,
		lon: -94.70365906334456,
		pop: 8087
	},
	{
		city: 'Fowler',
		state: 'Kansas',
		lat: 37.382248,
		lon: -100.195304782078864,
		pop: 590
	},
	{
		city: 'Frankfort',
		state: 'Kansas',
		lat: 39.703202,
		lon: -96.41700932944326,
		pop: 726
	},
	{
		city: 'Franklin',
		state: 'Kansas',
		lat: 37.523134,
		lon: -94.698807629120864,
		pop: 375
	},
	{
		city: 'Frederick',
		state: 'Kansas',
		lat: 38.5123975,
		lon: -98.267185014930334,
		pop: 18
	},
	{
		city: 'Fredonia',
		state: 'Kansas',
		lat: 37.5354935,
		lon: -95.825720365692547,
		pop: 2482
	},
	{
		city: 'Freeport',
		state: 'Kansas',
		lat: 37.1986255,
		lon: -97.854103133567335,
		pop: 5
	},
	{
		city: 'Frontenac',
		state: 'Kansas',
		lat: 37.4622785,
		lon: -94.687250601514648,
		pop: 3437
	},
	{
		city: 'Fulton',
		state: 'Kansas',
		lat: 38.009985,
		lon: -94.719522144431153,
		pop: 163
	},
	{
		city: 'Galatia',
		state: 'Kansas',
		lat: 38.6402935,
		lon: -98.957074025460827,
		pop: 39
	},
	{
		city: 'Galena',
		state: 'Kansas',
		lat: 37.0748235,
		lon: -94.636072787089859,
		pop: 3085
	},
	{
		city: 'Galesburg',
		state: 'Kansas',
		lat: 37.4719265,
		lon: -95.356526097122298,
		pop: 126
	},
	{
		city: 'Galva',
		state: 'Kansas',
		lat: 38.384344,
		lon: -97.537405659883746,
		pop: 870
	},
	{
		city: 'Garden City',
		state: 'Kansas',
		lat: 37.979228,
		lon: -100.864538557893212,
		pop: 26658
	},
	{
		city: 'Garden Plain',
		state: 'Kansas',
		lat: 37.6609205,
		lon: -97.682928500062616,
		pop: 849
	},
	{
		city: 'Gardner',
		state: 'Kansas',
		lat: 38.814507,
		lon: -94.929757020885688,
		pop: 19123
	},
	{
		city: 'Garfield',
		state: 'Kansas',
		lat: 38.075882,
		lon: -99.244088666980332,
		pop: 190
	},
	{
		city: 'Garnett',
		state: 'Kansas',
		lat: 38.286178,
		lon: -95.238524563817805,
		pop: 3415
	},
	{
		city: 'Gas',
		state: 'Kansas',
		lat: 37.9230775,
		lon: -95.344970773120281,
		pop: 564
	},
	{
		city: 'Gaylord',
		state: 'Kansas',
		lat: 39.6462515,
		lon: -98.846297406949347,
		pop: 114
	},
	{
		city: 'Gem',
		state: 'Kansas',
		lat: 39.4269365,
		lon: -100.896226626153833,
		pop: 88
	},
	{
		city: 'Geneseo',
		state: 'Kansas',
		lat: 38.5166915,
		lon: -98.155165947423399,
		pop: 267
	},
	{
		city: 'Geuda Springs',
		state: 'Kansas',
		lat: 37.110295,
		lon: -97.14851600028112,
		pop: 185
	},
	{
		city: 'Girard',
		state: 'Kansas',
		lat: 37.5104865,
		lon: -94.845326271243465,
		pop: 2789
	},
	{
		city: 'Glade',
		state: 'Kansas',
		lat: 39.682762,
		lon: -99.311242689010143,
		pop: 96
	},
	{
		city: 'Glasco',
		state: 'Kansas',
		lat: 39.360328,
		lon: -97.837402044035414,
		pop: 498
	},
	{
		city: 'Glen Elder',
		state: 'Kansas',
		lat: 39.5004685,
		lon: -98.306367090547937,
		pop: 445
	},
	{
		city: 'Goddard',
		state: 'Kansas',
		lat: 37.665195,
		lon: -97.562003409380281,
		pop: 4344
	},
	{
		city: 'Goessel',
		state: 'Kansas',
		lat: 38.2458705,
		lon: -97.346410465740291,
		pop: 539
	},
	{
		city: 'Goff',
		state: 'Kansas',
		lat: 39.6642055,
		lon: -95.93236802939856,
		pop: 126
	},
	{
		city: 'Goodland',
		state: 'Kansas',
		lat: 39.351212,
		lon: -101.711915614685864,
		pop: 4489
	},
	{
		city: 'Gorham',
		state: 'Kansas',
		lat: 38.8804885,
		lon: -99.02360897757984,
		pop: 334
	},
	{
		city: 'Gove City',
		state: 'Kansas',
		lat: 38.9592925,
		lon: -100.48708675,
		pop: 80
	},
	{
		city: 'Grainfield',
		state: 'Kansas',
		lat: 39.114803,
		lon: -100.468943290912677,
		pop: 277
	},
	{
		city: 'Grandview Plaza',
		state: 'Kansas',
		lat: 39.0328505,
		lon: -96.796410540253433,
		pop: 1560
	},
	{
		city: 'Grantville',
		state: 'Kansas',
		lat: 39.0803375,
		lon: -95.556988116276841,
		pop: 180
	},
	{
		city: 'Great Bend',
		state: 'Kansas',
		lat: 38.3602325,
		lon: -98.775137975363023,
		pop: 15995
	},
	{
		city: 'Greeley',
		state: 'Kansas',
		lat: 38.3680825,
		lon: -95.127313051063567,
		pop: 302
	},
	{
		city: 'Green',
		state: 'Kansas',
		lat: 39.4303345,
		lon: -97.000658693668527,
		pop: 128
	},
	{
		city: 'Greenleaf',
		state: 'Kansas',
		lat: 39.7271405,
		lon: -96.97928650855097,
		pop: 331
	},
	{
		city: 'Greensburg',
		state: 'Kansas',
		lat: 37.604834,
		lon: -99.292514969954368,
		pop: 777
	},
	{
		city: 'Grenola',
		state: 'Kansas',
		lat: 37.3496915,
		lon: -96.449064594784318,
		pop: 216
	},
	{
		city: 'Gridley',
		state: 'Kansas',
		lat: 38.0972615,
		lon: -95.883566000070459,
		pop: 341
	},
	{
		city: 'Grinnell',
		state: 'Kansas',
		lat: 39.125673,
		lon: -100.631877969424465,
		pop: 259
	},
	{
		city: 'Gypsum',
		state: 'Kansas',
		lat: 38.7056045,
		lon: -97.426501861856025,
		pop: 405
	},
	{
		city: 'Haddam',
		state: 'Kansas',
		lat: 39.8553335,
		lon: -97.305866438718652,
		pop: 104
	},
	{
		city: 'Halstead',
		state: 'Kansas',
		lat: 37.9987235,
		lon: -97.509494329957789,
		pop: 2085
	},
	{
		city: 'Hamilton',
		state: 'Kansas',
		lat: 37.9812865,
		lon: -96.164256824851691,
		pop: 268
	},
	{
		city: 'Hamlin',
		state: 'Kansas',
		lat: 39.9156245,
		lon: -95.627854562117278,
		pop: 46
	},
	{
		city: 'Hanover',
		state: 'Kansas',
		lat: 39.893147,
		lon: -96.877350095572353,
		pop: 682
	},
	{
		city: 'Hanston',
		state: 'Kansas',
		lat: 38.1215785,
		lon: -99.712308905206868,
		pop: 206
	},
	{
		city: 'Hardtner',
		state: 'Kansas',
		lat: 37.015406,
		lon: -98.649326261510382,
		pop: 172
	},
	{
		city: 'Harper',
		state: 'Kansas',
		lat: 37.283378,
		lon: -98.025345227495109,
		pop: 1473
	},
	{
		city: 'Harris',
		state: 'Kansas',
		lat: 38.322442,
		lon: -95.448589168788516,
		pop: 51
	},
	{
		city: 'Hartford',
		state: 'Kansas',
		lat: 38.3085315,
		lon: -95.956275348117899,
		pop: 371
	},
	{
		city: 'Harveyville',
		state: 'Kansas',
		lat: 38.789248,
		lon: -95.962448660800447,
		pop: 236
	},
	{
		city: 'Havana',
		state: 'Kansas',
		lat: 37.0923595,
		lon: -95.942833314131747,
		pop: 104
	},
	{
		city: 'Haven',
		state: 'Kansas',
		lat: 37.901843,
		lon: -97.782020256506513,
		pop: 1237
	},
	{
		city: 'Havensville',
		state: 'Kansas',
		lat: 39.511668,
		lon: -96.076234644990819,
		pop: 133
	},
	{
		city: 'Haviland',
		state: 'Kansas',
		lat: 37.6170495,
		lon: -99.10612877761011,
		pop: 701
	},
	{
		city: 'Hays',
		state: 'Kansas',
		lat: 38.885016,
		lon: -99.320131641613926,
		pop: 20510
	},
	{
		city: 'Haysville',
		state: 'Kansas',
		lat: 37.5653255,
		lon: -97.352747768658645,
		pop: 10826
	},
	{
		city: 'Hazelton',
		state: 'Kansas',
		lat: 37.0890335,
		lon: -98.402809318782303,
		pop: 93
	},
	{
		city: 'Healy',
		state: 'Kansas',
		lat: 38.6026535,
		lon: -100.613820143564254,
		pop: 234
	},
	{
		city: 'Hepler',
		state: 'Kansas',
		lat: 37.6625935,
		lon: -94.969632363324081,
		pop: 132
	},
	{
		city: 'Herington',
		state: 'Kansas',
		lat: 38.6961665,
		lon: -96.797828146690136,
		pop: 2526
	},
	{
		city: 'Herndon',
		state: 'Kansas',
		lat: 39.908255,
		lon: -100.786198608696537,
		pop: 129
	},
	{
		city: 'Hesston',
		state: 'Kansas',
		lat: 38.1406945,
		lon: -97.417874546773078,
		pop: 3709
	},
	{
		city: 'Hiawatha',
		state: 'Kansas',
		lat: 39.8536135,
		lon: -95.536184808888621,
		pop: 3172
	},
	{
		city: 'Highland',
		state: 'Kansas',
		lat: 39.8600185,
		lon: -95.265089766931212,
		pop: 1012
	},
	{
		city: 'Hill City',
		state: 'Kansas',
		lat: 39.367241,
		lon: -99.847658376764628,
		pop: 1474
	},
	{
		city: 'Hillsboro',
		state: 'Kansas',
		lat: 38.348608,
		lon: -97.201027352390838,
		pop: 2993
	},
	{
		city: 'Hillsdale',
		state: 'Kansas',
		lat: 38.6555175,
		lon: -94.85875752062482,
		pop: 229
	},
	{
		city: 'Hoisington',
		state: 'Kansas',
		lat: 38.518943,
		lon: -98.777373487903873,
		pop: 2706
	},
	{
		city: 'Holcomb',
		state: 'Kansas',
		lat: 37.990097,
		lon: -100.98982790013774,
		pop: 2094
	},
	{
		city: 'Hollenberg',
		state: 'Kansas',
		lat: 39.980902,
		lon: -96.992413755249345,
		pop: 21
	},
	{
		city: 'Holton',
		state: 'Kansas',
		lat: 39.464272,
		lon: -95.741672286259259,
		pop: 3329
	},
	{
		city: 'Holyrood',
		state: 'Kansas',
		lat: 38.5884905,
		lon: -98.411239373472938,
		pop: 447
	},
	{
		city: 'Home',
		state: 'Kansas',
		lat: 39.8414905,
		lon: -96.518701483455885,
		pop: 160
	},
	{
		city: 'Hope',
		state: 'Kansas',
		lat: 38.690979,
		lon: -97.078241090870307,
		pop: 368
	},
	{
		city: 'Horace',
		state: 'Kansas',
		lat: 38.477075,
		lon: -101.790344336709211,
		pop: 70
	},
	{
		city: 'Horton',
		state: 'Kansas',
		lat: 39.6656365,
		lon: -95.529351231022815,
		pop: 1776
	},
	{
		city: 'Howard',
		state: 'Kansas',
		lat: 37.469347,
		lon: -96.262314403846148,
		pop: 687
	},
	{
		city: 'Hoxie',
		state: 'Kansas',
		lat: 39.355128,
		lon: -100.440246422381307,
		pop: 1201
	},
	{
		city: 'Hoyt',
		state: 'Kansas',
		lat: 39.2497295,
		lon: -95.705987437484964,
		pop: 669
	},
	{
		city: 'Hudson',
		state: 'Kansas',
		lat: 38.103438,
		lon: -98.660156591926182,
		pop: 129
	},
	{
		city: 'Hugoton',
		state: 'Kansas',
		lat: 37.175422,
		lon: -101.347669594955619,
		pop: 3904
	},
	{
		city: 'Humboldt',
		state: 'Kansas',
		lat: 37.812073,
		lon: -95.437950799078322,
		pop: 1953
	},
	{
		city: 'Hunnewell',
		state: 'Kansas',
		lat: 37.0047635,
		lon: -97.407712484244996,
		pop: 67
	},
	{
		city: 'Hunter',
		state: 'Kansas',
		lat: 39.234691,
		lon: -98.395034454268284,
		pop: 57
	},
	{
		city: 'Huron',
		state: 'Kansas',
		lat: 39.643909,
		lon: -95.346968242353384,
		pop: 54
	},
	{
		city: 'Hutchinson',
		state: 'Kansas',
		lat: 38.063636,
		lon: -97.925660699972241,
		pop: 42080
	},
	{
		city: 'Independence',
		state: 'Kansas',
		lat: 37.155634,
		lon: -95.781554301355015,
		pop: 9483
	},
	{
		city: 'Ingalls',
		state: 'Kansas',
		lat: 37.828252,
		lon: -100.452562600706528,
		pop: 306
	},
	{
		city: 'Inman',
		state: 'Kansas',
		lat: 38.230954,
		lon: -97.770958578042624,
		pop: 1377
	},
	{
		city: 'Iola',
		state: 'Kansas',
		lat: 37.9256905,
		lon: -95.395379911983468,
		pop: 5704
	},
	{
		city: 'Isabel',
		state: 'Kansas',
		lat: 37.4672735,
		lon: -98.551510676288657,
		pop: 90
	},
	{
		city: 'Iuka',
		state: 'Kansas',
		lat: 37.729766,
		lon: -98.730124700554512,
		pop: 163
	},
	{
		city: 'Jamestown',
		state: 'Kansas',
		lat: 39.599051,
		lon: -97.861178632944217,
		pop: 286
	},
	{
		city: 'Jennings',
		state: 'Kansas',
		lat: 39.6816935,
		lon: -100.293583359225039,
		pop: 96
	},
	{
		city: 'Jetmore',
		state: 'Kansas',
		lat: 38.0826535,
		lon: -99.894178745080865,
		pop: 867
	},
	{
		city: 'Jewell',
		state: 'Kansas',
		lat: 39.6719865,
		lon: -98.1519761631356,
		pop: 432
	},
	{
		city: 'Johnson City',
		state: 'Kansas',
		lat: 37.56891,
		lon: -101.752691364897757,
		pop: 1495
	},
	{
		city: 'Junction City',
		state: 'Kansas',
		lat: 39.025818,
		lon: -96.839052011517964,
		pop: 23353
	},
	{
		city: 'Kanopolis',
		state: 'Kansas',
		lat: 38.707434,
		lon: -98.158182177476249,
		pop: 492
	},
	{
		city: 'Kanorado',
		state: 'Kansas',
		lat: 39.336993,
		lon: -102.037414118441461,
		pop: 153
	},
	{
		city: 'Kansas City',
		state: 'Kansas',
		lat: 39.1230125,
		lon: -94.756838130084745,
		pop: 145786
	},
	{
		city: 'Kechi',
		state: 'Kansas',
		lat: 37.800905,
		lon: -97.280889299023414,
		pop: 1909
	},
	{
		city: 'Kensington',
		state: 'Kansas',
		lat: 39.76706,
		lon: -99.033061927210142,
		pop: 473
	},
	{
		city: 'Kickapoo Site 1',
		state: 'Kansas',
		lat: 39.715169,
		lon: -95.648402657442233,
		pop: 101
	},
	{
		city: 'Kickapoo Site 2',
		state: 'Kansas',
		lat: 39.703724,
		lon: -95.650332387018011,
		pop: 34
	},
	{
		city: 'Kickapoo Site 5',
		state: 'Kansas',
		lat: 39.6746825,
		lon: -95.681585597014873,
		pop: 66
	},
	{
		city: 'Kickapoo Site 6',
		state: 'Kansas',
		lat: 39.690548,
		lon: -95.690400605425395,
		pop: 15
	},
	{
		city: 'Kickapoo Site 7',
		state: 'Kansas',
		lat: 39.6889135,
		lon: -95.672197177996566,
		pop: 66
	},
	{
		city: 'Kickapoo Tribal Center',
		state: 'Kansas',
		lat: 39.6708515,
		lon: -95.645290574116046,
		pop: 194
	},
	{
		city: 'Kincaid',
		state: 'Kansas',
		lat: 38.0811985,
		lon: -95.155245760638294,
		pop: 122
	},
	{
		city: 'Kingman',
		state: 'Kansas',
		lat: 37.648145,
		lon: -98.108559678809598,
		pop: 3177
	},
	{
		city: 'Kinsley',
		state: 'Kansas',
		lat: 37.923634,
		lon: -99.410264500077432,
		pop: 1457
	},
	{
		city: 'Kiowa',
		state: 'Kansas',
		lat: 37.018567,
		lon: -98.485278425355432,
		pop: 1026
	},
	{
		city: 'Kipp',
		state: 'Kansas',
		lat: 38.7843,
		lon: -97.455317180882346,
		pop: 59
	},
	{
		city: 'Kirwin',
		state: 'Kansas',
		lat: 39.6685505,
		lon: -99.122458801300482,
		pop: 171
	},
	{
		city: 'Kismet',
		state: 'Kansas',
		lat: 37.2042725,
		lon: -100.702763661501024,
		pop: 459
	},
	{
		city: 'La Crosse',
		state: 'Kansas',
		lat: 38.533184,
		lon: -99.30964975560893,
		pop: 1342
	},
	{
		city: 'La Cygne',
		state: 'Kansas',
		lat: 38.348753,
		lon: -94.756586143143537,
		pop: 1149
	},
	{
		city: 'La Harpe',
		state: 'Kansas',
		lat: 37.915536,
		lon: -95.303318327951985,
		pop: 578
	},
	{
		city: 'Labette',
		state: 'Kansas',
		lat: 37.230552,
		lon: -95.18190492231949,
		pop: 78
	},
	{
		city: 'Lake Quivira',
		state: 'Kansas',
		lat: 39.039319,
		lon: -94.770169561541337,
		pop: 906
	},
	{
		city: 'Lakin',
		state: 'Kansas',
		lat: 37.9393815,
		lon: -101.259608461538463,
		pop: 2216
	},
	{
		city: 'Lancaster',
		state: 'Kansas',
		lat: 39.571416,
		lon: -95.304316567098908,
		pop: 298
	},
	{
		city: 'Lane',
		state: 'Kansas',
		lat: 38.438904,
		lon: -95.083936728908185,
		pop: 225
	},
	{
		city: 'Langdon',
		state: 'Kansas',
		lat: 37.8526315,
		lon: -98.324815996738366,
		pop: 42
	},
	{
		city: 'Lansing',
		state: 'Kansas',
		lat: 39.2396135,
		lon: -94.900328372567714,
		pop: 11265
	},
	{
		city: 'Larned',
		state: 'Kansas',
		lat: 38.1831305,
		lon: -99.10251759789071,
		pop: 4054
	},
	{
		city: 'Latham',
		state: 'Kansas',
		lat: 37.5350875,
		lon: -96.642689469738798,
		pop: 139
	},
	{
		city: 'Latimer',
		state: 'Kansas',
		lat: 38.7386435,
		lon: -96.845574048779696,
		pop: 20
	},
	{
		city: 'Lawrence',
		state: 'Kansas',
		lat: 38.9736425,
		lon: -95.273492487141141,
		pop: 87643
	},
	{
		city: 'LeRoy',
		state: 'Kansas',
		lat: 38.085944,
		lon: -95.634372611895913,
		pop: 561
	},
	{
		city: 'Leavenworth',
		state: 'Kansas',
		lat: 39.325658,
		lon: -94.932462602428387,
		pop: 35251
	},
	{
		city: 'Leawood',
		state: 'Kansas',
		lat: 38.915013,
		lon: -94.628525805847417,
		pop: 31867
	},
	{
		city: 'Lebanon',
		state: 'Kansas',
		lat: 39.811094,
		lon: -98.557288654026962,
		pop: 218
	},
	{
		city: 'Lebo',
		state: 'Kansas',
		lat: 38.4155325,
		lon: -95.857942741943191,
		pop: 940
	},
	{
		city: 'Lecompton',
		state: 'Kansas',
		lat: 39.0353725,
		lon: -95.388052170885771,
		pop: 625
	},
	{
		city: 'Lehigh',
		state: 'Kansas',
		lat: 38.374916,
		lon: -97.302887657320866,
		pop: 175
	},
	{
		city: 'Lenexa',
		state: 'Kansas',
		lat: 38.9543695,
		lon: -94.816152033312505,
		pop: 48190
	},
	{
		city: 'Lenora',
		state: 'Kansas',
		lat: 39.61083,
		lon: -100.002783906777054,
		pop: 250
	},
	{
		city: 'Leon',
		state: 'Kansas',
		lat: 37.690518,
		lon: -96.781975716666238,
		pop: 704
	},
	{
		city: 'Leona',
		state: 'Kansas',
		lat: 39.78567,
		lon: -95.321378777119321,
		pop: 48
	},
	{
		city: 'Leonardville',
		state: 'Kansas',
		lat: 39.364106,
		lon: -96.858780470588229,
		pop: 449
	},
	{
		city: 'Leoti',
		state: 'Kansas',
		lat: 38.4847535,
		lon: -101.357552750496239,
		pop: 1534
	},
	{
		city: 'Levant',
		state: 'Kansas',
		lat: 39.388046,
		lon: -101.194726405208655,
		pop: 61
	},
	{
		city: 'Lewis',
		state: 'Kansas',
		lat: 37.937316,
		lon: -99.254618131500735,
		pop: 451
	},
	{
		city: 'Liberal',
		state: 'Kansas',
		lat: 37.0473335,
		lon: -100.919199371405966,
		pop: 20525
	},
	{
		city: 'Liberty',
		state: 'Kansas',
		lat: 37.15739,
		lon: -95.598401335897421,
		pop: 123
	},
	{
		city: 'Liebenthal',
		state: 'Kansas',
		lat: 38.653709,
		lon: -99.31792788168498,
		pop: 103
	},
	{
		city: 'Lincoln Center',
		state: 'Kansas',
		lat: 39.044098,
		lon: -98.146250136051179,
		pop: 1297
	},
	{
		city: 'Lincolnville',
		state: 'Kansas',
		lat: 38.4937435,
		lon: -96.961755782705467,
		pop: 203
	},
	{
		city: 'Lindsborg',
		state: 'Kansas',
		lat: 38.5772605,
		lon: -97.674990868703361,
		pop: 3458
	},
	{
		city: 'Linn',
		state: 'Kansas',
		lat: 39.679752,
		lon: -97.086320690062621,
		pop: 804
	},
	{
		city: 'Linn Valley',
		state: 'Kansas',
		lat: 38.375283,
		lon: -94.714119751759327,
		pop: 410
	},
	{
		city: 'Linwood',
		state: 'Kansas',
		lat: 39.0032385,
		lon: -95.036822782962716,
		pop: 375
	},
	{
		city: 'Little River',
		state: 'Kansas',
		lat: 38.396872,
		lon: -98.013228215581037,
		pop: 557
	},
	{
		city: 'Logan',
		state: 'Kansas',
		lat: 39.661562,
		lon: -99.56711524512987,
		pop: 589
	},
	{
		city: 'Lone Elm',
		state: 'Kansas',
		lat: 38.0797755,
		lon: -95.243039495779641,
		pop: 25
	},
	{
		city: 'Long Island',
		state: 'Kansas',
		lat: 39.9460815,
		lon: -99.535973979012141,
		pop: 134
	},
	{
		city: 'Longford',
		state: 'Kansas',
		lat: 39.17283,
		lon: -97.32860261958804,
		pop: 79
	},
	{
		city: 'Longton',
		state: 'Kansas',
		lat: 37.377586,
		lon: -96.080344506142964,
		pop: 348
	},
	{
		city: 'Lorraine',
		state: 'Kansas',
		lat: 38.569266,
		lon: -98.317125400960506,
		pop: 138
	},
	{
		city: 'Lost Springs',
		state: 'Kansas',
		lat: 38.5674945,
		lon: -96.965850143052194,
		pop: 70
	},
	{
		city: 'Louisburg',
		state: 'Kansas',
		lat: 38.6209995,
		lon: -94.687480833269348,
		pop: 4315
	},
	{
		city: 'Louisville',
		state: 'Kansas',
		lat: 39.249552,
		lon: -96.314508740571426,
		pop: 188
	},
	{
		city: 'Lowell',
		state: 'Kansas',
		lat: 37.051188,
		lon: -94.703567316961426,
		pop: 283
	},
	{
		city: 'Lucas',
		state: 'Kansas',
		lat: 39.0594715,
		lon: -98.540686454310332,
		pop: 393
	},
	{
		city: 'Luray',
		state: 'Kansas',
		lat: 39.1146565,
		lon: -98.69182834708738,
		pop: 194
	},
	{
		city: 'Lyndon',
		state: 'Kansas',
		lat: 38.6126795,
		lon: -95.685470548401639,
		pop: 1052
	},
	{
		city: 'Lyons',
		state: 'Kansas',
		lat: 38.3469505,
		lon: -98.20765856062053,
		pop: 3739
	},
	{
		city: 'Macksville',
		state: 'Kansas',
		lat: 37.959241,
		lon: -98.969181943304093,
		pop: 549
	},
	{
		city: 'Madison',
		state: 'Kansas',
		lat: 38.1325245,
		lon: -96.135932386133788,
		pop: 701
	},
	{
		city: 'Mahaska',
		state: 'Kansas',
		lat: 39.9880095,
		lon: -97.353465381481357,
		pop: 83
	},
	{
		city: 'Maize',
		state: 'Kansas',
		lat: 37.766071,
		lon: -97.466041889588865,
		pop: 3420
	},
	{
		city: 'Manchester',
		state: 'Kansas',
		lat: 39.0906535,
		lon: -97.321066327987552,
		pop: 95
	},
	{
		city: 'Manhattan',
		state: 'Kansas',
		lat: 39.1910025,
		lon: -96.59236755453307,
		pop: 52281
	},
	{
		city: 'Mankato',
		state: 'Kansas',
		lat: 39.786784,
		lon: -98.205433206977716,
		pop: 869
	},
	{
		city: 'Manter',
		state: 'Kansas',
		lat: 37.523418,
		lon: -101.884469101967198,
		pop: 171
	},
	{
		city: 'Maple Hill',
		state: 'Kansas',
		lat: 39.0862195,
		lon: -96.028133426464066,
		pop: 620
	},
	{
		city: 'Mapleton',
		state: 'Kansas',
		lat: 38.0154665,
		lon: -94.883418984123438,
		pop: 84
	},
	{
		city: 'Marienthal',
		state: 'Kansas',
		lat: 38.4890235,
		lon: -101.220018780465324,
		pop: 71
	},
	{
		city: 'Marion',
		state: 'Kansas',
		lat: 38.3383735,
		lon: -96.992057306730658,
		pop: 1927
	},
	{
		city: 'Marquette',
		state: 'Kansas',
		lat: 38.557191,
		lon: -97.834562456547758,
		pop: 641
	},
	{
		city: 'Marysville',
		state: 'Kansas',
		lat: 39.84632,
		lon: -96.64189166782927,
		pop: 3294
	},
	{
		city: 'Matfield Green',
		state: 'Kansas',
		lat: 38.15862,
		lon: -96.562161064658483,
		pop: 47
	},
	{
		city: 'Mayetta',
		state: 'Kansas',
		lat: 39.338221,
		lon: -95.72125420385305,
		pop: 341
	},
	{
		city: 'Mayfield',
		state: 'Kansas',
		lat: 37.261794,
		lon: -97.546449642992428,
		pop: 113
	},
	{
		city: 'McConnell AFB',
		state: 'Kansas',
		lat: 37.618995,
		lon: -97.26074975848266,
		pop: 1777
	},
	{
		city: 'McCracken',
		state: 'Kansas',
		lat: 38.5820295,
		lon: -99.568991299189804,
		pop: 190
	},
	{
		city: 'McCune',
		state: 'Kansas',
		lat: 37.353109,
		lon: -95.019429484126988,
		pop: 405
	},
	{
		city: 'McDonald',
		state: 'Kansas',
		lat: 39.7848965,
		lon: -101.369365450199211,
		pop: 160
	},
	{
		city: 'McFarland',
		state: 'Kansas',
		lat: 39.0543815,
		lon: -96.237441431799596,
		pop: 256
	},
	{
		city: 'McLouth',
		state: 'Kansas',
		lat: 39.1963245,
		lon: -95.207200351154611,
		pop: 880
	},
	{
		city: 'McPherson',
		state: 'Kansas',
		lat: 38.370327,
		lon: -97.653542314793583,
		pop: 13155
	},
	{
		city: 'Meade',
		state: 'Kansas',
		lat: 37.2842595,
		lon: -100.3403083198166,
		pop: 1721
	},
	{
		city: 'Medicine Lodge',
		state: 'Kansas',
		lat: 37.285633,
		lon: -98.584347360978924,
		pop: 2009
	},
	{
		city: 'Melvern',
		state: 'Kansas',
		lat: 38.506588,
		lon: -95.636842918467579,
		pop: 385
	},
	{
		city: 'Menlo',
		state: 'Kansas',
		lat: 39.35595,
		lon: -100.724408579596414,
		pop: 61
	},
	{
		city: 'Meriden',
		state: 'Kansas',
		lat: 39.189751,
		lon: -95.568842138448545,
		pop: 813
	},
	{
		city: 'Merriam',
		state: 'Kansas',
		lat: 39.018079,
		lon: -94.690997396122697,
		pop: 11003
	},
	{
		city: 'Milan',
		state: 'Kansas',
		lat: 37.258506,
		lon: -97.673938090039428,
		pop: 82
	},
	{
		city: 'Mildred',
		state: 'Kansas',
		lat: 38.025742,
		lon: -95.171127738039559,
		pop: 28
	},
	{
		city: 'Milford',
		state: 'Kansas',
		lat: 39.1715395,
		lon: -96.911768360938254,
		pop: 530
	},
	{
		city: 'Milton',
		state: 'Kansas',
		lat: 37.441586,
		lon: -97.769582146111631,
		pop: 155
	},
	{
		city: 'Miltonvale',
		state: 'Kansas',
		lat: 39.3498905,
		lon: -97.451402205606854,
		pop: 539
	},
	{
		city: 'Minneapolis',
		state: 'Kansas',
		lat: 39.119867,
		lon: -97.703386196394121,
		pop: 2032
	},
	{
		city: 'Minneola',
		state: 'Kansas',
		lat: 37.4420315,
		lon: -100.013702571163236,
		pop: 745
	},
	{
		city: 'Mission',
		state: 'Kansas',
		lat: 39.026499,
		lon: -94.658655026090614,
		pop: 3498
	},
	{
		city: 'Mission Hills',
		state: 'Kansas',
		lat: 39.0170875,
		lon: -94.616156963897339,
		pop: 178
	},
	{
		city: 'Mission Woods',
		state: 'Kansas',
		lat: 39.0336435,
		lon: -94.610050145030414,
		pop: 9323
	},
	{
		city: 'Moline',
		state: 'Kansas',
		lat: 37.364653,
		lon: -96.302764933925801,
		pop: 371
	},
	{
		city: 'Montezuma',
		state: 'Kansas',
		lat: 37.5961225,
		lon: -100.443471960526324,
		pop: 966
	},
	{
		city: 'Moran',
		state: 'Kansas',
		lat: 37.9168335,
		lon: -95.172245134884278,
		pop: 558
	},
	{
		city: 'Morganville',
		state: 'Kansas',
		lat: 39.4665915,
		lon: -97.2044565269766,
		pop: 192
	},
	{
		city: 'Morland',
		state: 'Kansas',
		lat: 39.349652,
		lon: -100.075861797049441,
		pop: 154
	},
	{
		city: 'Morrill',
		state: 'Kansas',
		lat: 39.9294415,
		lon: -95.694621662499998,
		pop: 230
	},
	{
		city: 'Morrowville',
		state: 'Kansas',
		lat: 39.844034,
		lon: -97.172945151150046,
		pop: 155
	},
	{
		city: 'Moscow',
		state: 'Kansas',
		lat: 37.3273355,
		lon: -101.205571383630001,
		pop: 310
	},
	{
		city: 'Mound City',
		state: 'Kansas',
		lat: 38.1431065,
		lon: -94.811889574123867,
		pop: 694
	},
	{
		city: 'Mound Valley',
		state: 'Kansas',
		lat: 37.2063225,
		lon: -95.403309712133279,
		pop: 407
	},
	{
		city: 'Moundridge',
		state: 'Kansas',
		lat: 38.2056975,
		lon: -97.506184010187241,
		pop: 1737
	},
	{
		city: 'Mount Hope',
		state: 'Kansas',
		lat: 37.872347,
		lon: -97.65792340797546,
		pop: 813
	},
	{
		city: 'Mulberry',
		state: 'Kansas',
		lat: 37.5552115,
		lon: -94.626161991637474,
		pop: 520
	},
	{
		city: 'Mullinville',
		state: 'Kansas',
		lat: 37.588114,
		lon: -99.475763497836837,
		pop: 255
	},
	{
		city: 'Mulvane',
		state: 'Kansas',
		lat: 37.4816815,
		lon: -97.238261325555555,
		pop: 6111
	},
	{
		city: 'Munden',
		state: 'Kansas',
		lat: 39.911571,
		lon: -97.538405448929581,
		pop: 100
	},
	{
		city: 'Munjor',
		state: 'Kansas',
		lat: 38.8121775,
		lon: -99.266810177075541,
		pop: 213
	},
	{
		city: 'Muscotah',
		state: 'Kansas',
		lat: 39.5530275,
		lon: -95.519998435376991,
		pop: 176
	},
	{
		city: 'Narka',
		state: 'Kansas',
		lat: 39.959725,
		lon: -97.426809850588,
		pop: 94
	},
	{
		city: 'Nashville',
		state: 'Kansas',
		lat: 37.438802,
		lon: -98.421406100875402,
		pop: 64
	},
	{
		city: 'Natoma',
		state: 'Kansas',
		lat: 39.186901,
		lon: -99.024216618535178,
		pop: 335
	},
	{
		city: 'Neodesha',
		state: 'Kansas',
		lat: 37.4104605,
		lon: -95.67060750987082,
		pop: 2486
	},
	{
		city: 'Neosho Falls',
		state: 'Kansas',
		lat: 38.0057915,
		lon: -95.555693539560352,
		pop: 141
	},
	{
		city: 'Neosho Rapids',
		state: 'Kansas',
		lat: 38.368423,
		lon: -95.99197445276431,
		pop: 265
	},
	{
		city: 'Ness City',
		state: 'Kansas',
		lat: 38.4533495,
		lon: -99.90508648616769,
		pop: 1449
	},
	{
		city: 'Netawaka',
		state: 'Kansas',
		lat: 39.602818,
		lon: -95.720522677308054,
		pop: 143
	},
	{
		city: 'New Albany',
		state: 'Kansas',
		lat: 37.5667155,
		lon: -95.931757617963143,
		pop: 56
	},
	{
		city: 'New Cambria',
		state: 'Kansas',
		lat: 38.8788185,
		lon: -97.506737092776035,
		pop: 126
	},
	{
		city: 'New Strawn',
		state: 'Kansas',
		lat: 38.261424,
		lon: -95.741799726187196,
		pop: 394
	},
	{
		city: 'Newton',
		state: 'Kansas',
		lat: 38.0372665,
		lon: -97.349822281288283,
		pop: 19132
	},
	{
		city: 'Nickerson',
		state: 'Kansas',
		lat: 38.146135,
		lon: -98.088616970628408,
		pop: 1070
	},
	{
		city: 'Niotaze',
		state: 'Kansas',
		lat: 37.0684975,
		lon: -96.014096273063416,
		pop: 82
	},
	{
		city: 'Norcatur',
		state: 'Kansas',
		lat: 39.8348305,
		lon: -100.188764940704473,
		pop: 151
	},
	{
		city: 'North Newton',
		state: 'Kansas',
		lat: 38.075993,
		lon: -97.345836673038221,
		pop: 1759
	},
	{
		city: 'Norton',
		state: 'Kansas',
		lat: 39.8381815,
		lon: -99.893480918520098,
		pop: 2928
	},
	{
		city: 'Nortonville',
		state: 'Kansas',
		lat: 39.414161,
		lon: -95.329948123156299,
		pop: 637
	},
	{
		city: 'Norwich',
		state: 'Kansas',
		lat: 37.456789,
		lon: -97.847616235624031,
		pop: 491
	},
	{
		city: 'Oak Hill',
		state: 'Kansas',
		lat: 39.2465635,
		lon: -97.342993480749215,
		pop: 24
	},
	{
		city: 'Oaklawn-Sunview',
		state: 'Kansas',
		lat: 37.608007,
		lon: -97.299883425847455,
		pop: 3276
	},
	{
		city: 'Oakley',
		state: 'Kansas',
		lat: 39.1230085,
		lon: -100.851852210669634,
		pop: 2045
	},
	{
		city: 'Oberlin',
		state: 'Kansas',
		lat: 39.824891,
		lon: -100.529661879742804,
		pop: 1788
	},
	{
		city: 'Odin',
		state: 'Kansas',
		lat: 38.565366,
		lon: -98.608799104235558,
		pop: 101
	},
	{
		city: 'Offerle',
		state: 'Kansas',
		lat: 37.8901135,
		lon: -99.560792572489788,
		pop: 199
	},
	{
		city: 'Ogden',
		state: 'Kansas',
		lat: 39.114818,
		lon: -96.710311318923317,
		pop: 2087
	},
	{
		city: 'Oketo',
		state: 'Kansas',
		lat: 39.9630025,
		lon: -96.599088638186572,
		pop: 66
	},
	{
		city: 'Olathe',
		state: 'Kansas',
		lat: 38.8818835,
		lon: -94.804645292347061,
		pop: 125872
	},
	{
		city: 'Olivet',
		state: 'Kansas',
		lat: 38.4820865,
		lon: -95.752064624719793,
		pop: 67
	},
	{
		city: 'Olmitz',
		state: 'Kansas',
		lat: 38.5170495,
		lon: -98.936312369117104,
		pop: 114
	},
	{
		city: 'Olpe',
		state: 'Kansas',
		lat: 38.2623255,
		lon: -96.166556736784145,
		pop: 546
	},
	{
		city: 'Olsburg',
		state: 'Kansas',
		lat: 39.4313595,
		lon: -96.615133931245893,
		pop: 219
	},
	{
		city: 'Onaga',
		state: 'Kansas',
		lat: 39.488963,
		lon: -96.17008334378059,
		pop: 702
	},
	{
		city: 'Oneida',
		state: 'Kansas',
		lat: 39.866885,
		lon: -95.939782024283318,
		pop: 75
	},
	{
		city: 'Osage City',
		state: 'Kansas',
		lat: 38.629947,
		lon: -95.829110243980736,
		pop: 2943
	},
	{
		city: 'Osawatomie',
		state: 'Kansas',
		lat: 38.498679,
		lon: -94.955984974253084,
		pop: 4447
	},
	{
		city: 'Osborne',
		state: 'Kansas',
		lat: 39.4440875,
		lon: -98.699885274089695,
		pop: 1431
	},
	{
		city: 'Oskaloosa',
		state: 'Kansas',
		lat: 39.2151085,
		lon: -95.316887652828626,
		pop: 1113
	},
	{
		city: 'Oswego',
		state: 'Kansas',
		lat: 37.1685585,
		lon: -95.112747672968382,
		pop: 1829
	},
	{
		city: 'Otis',
		state: 'Kansas',
		lat: 38.534597,
		lon: -99.052808324687149,
		pop: 282
	},
	{
		city: 'Ottawa',
		state: 'Kansas',
		lat: 38.6076095,
		lon: -95.267222,
		pop: 12649
	},
	{
		city: 'Overbrook',
		state: 'Kansas',
		lat: 38.7807905,
		lon: -95.55652738851542,
		pop: 1058
	},
	{
		city: 'Overland Park',
		state: 'Kansas',
		lat: 38.902323,
		lon: -94.694201876548675,
		pop: 173372
	},
	{
		city: 'Oxford',
		state: 'Kansas',
		lat: 37.273495,
		lon: -97.170521928478365,
		pop: 1049
	},
	{
		city: 'Ozawkie',
		state: 'Kansas',
		lat: 39.235029,
		lon: -95.465465224332007,
		pop: 645
	},
	{
		city: 'Palco',
		state: 'Kansas',
		lat: 39.2544975,
		lon: -99.563653517557981,
		pop: 277
	},
	{
		city: 'Palmer',
		state: 'Kansas',
		lat: 39.633046,
		lon: -97.138543322321425,
		pop: 111
	},
	{
		city: 'Paola',
		state: 'Kansas',
		lat: 38.5795825,
		lon: -94.868706892511454,
		pop: 5602
	},
	{
		city: 'Paradise',
		state: 'Kansas',
		lat: 39.1141845,
		lon: -98.918447445615442,
		pop: 49
	},
	{
		city: 'Park',
		state: 'Kansas',
		lat: 39.111929,
		lon: -100.362690806525322,
		pop: 126
	},
	{
		city: 'Park City',
		state: 'Kansas',
		lat: 37.818855,
		lon: -97.321846761176715,
		pop: 7297
	},
	{
		city: 'Parker',
		state: 'Kansas',
		lat: 38.328124,
		lon: -94.990473554613544,
		pop: 277
	},
	{
		city: 'Parkerfield',
		state: 'Kansas',
		lat: 37.0661265,
		lon: -96.996240251712322,
		pop: 426
	},
	{
		city: 'Parkerville',
		state: 'Kansas',
		lat: 38.763794,
		lon: -96.661718816671566,
		pop: 59
	},
	{
		city: 'Parsons',
		state: 'Kansas',
		lat: 37.3473495,
		lon: -95.272956838227486,
		pop: 10500
	},
	{
		city: 'Partridge',
		state: 'Kansas',
		lat: 37.966072,
		lon: -98.092311335732973,
		pop: 248
	},
	{
		city: 'Pawnee Rock',
		state: 'Kansas',
		lat: 38.264839,
		lon: -98.98312200246562,
		pop: 252
	},
	{
		city: 'Paxico',
		state: 'Kansas',
		lat: 39.0691325,
		lon: -96.166679218331296,
		pop: 221
	},
	{
		city: 'Peabody',
		state: 'Kansas',
		lat: 38.168737,
		lon: -97.107059195532969,
		pop: 1210
	},
	{
		city: 'Penalosa',
		state: 'Kansas',
		lat: 37.7153825,
		lon: -98.318610010847792,
		pop: 17
	},
	{
		city: 'Perry',
		state: 'Kansas',
		lat: 39.0713355,
		lon: -95.383711382005757,
		pop: 929
	},
	{
		city: 'Peru',
		state: 'Kansas',
		lat: 37.0814955,
		lon: -96.09801725072974,
		pop: 139
	},
	{
		city: 'Phillipsburg',
		state: 'Kansas',
		lat: 39.755425,
		lon: -99.321947895926229,
		pop: 2581
	},
	{
		city: 'Piqua',
		state: 'Kansas',
		lat: 37.9220075,
		lon: -95.537496504822514,
		pop: 107
	},
	{
		city: 'Pittsburg',
		state: 'Kansas',
		lat: 37.416984,
		lon: -94.695715833373228,
		pop: 20233
	},
	{
		city: 'Plains',
		state: 'Kansas',
		lat: 37.26464,
		lon: -100.589763323931237,
		pop: 1146
	},
	{
		city: 'Plainville',
		state: 'Kansas',
		lat: 39.2335675,
		lon: -99.303171589467183,
		pop: 1903
	},
	{
		city: 'Pleasanton',
		state: 'Kansas',
		lat: 38.1908235,
		lon: -94.690977781154601,
		pop: 1216
	},
	{
		city: 'Plevna',
		state: 'Kansas',
		lat: 37.9713025,
		lon: -98.3093292288697,
		pop: 98
	},
	{
		city: 'Pomona',
		state: 'Kansas',
		lat: 38.611463,
		lon: -95.453050368010224,
		pop: 832
	},
	{
		city: 'Portis',
		state: 'Kansas',
		lat: 39.5632615,
		lon: -98.691324666975447,
		pop: 103
	},
	{
		city: 'Potwin',
		state: 'Kansas',
		lat: 37.938129,
		lon: -97.01768817538462,
		pop: 449
	},
	{
		city: 'Powhattan',
		state: 'Kansas',
		lat: 39.761156,
		lon: -95.63324791,
		pop: 77
	},
	{
		city: 'Prairie View',
		state: 'Kansas',
		lat: 39.8320405,
		lon: -99.573404740893977,
		pop: 134
	},
	{
		city: 'Prairie Village',
		state: 'Kansas',
		lat: 38.986624,
		lon: -94.633213560218977,
		pop: 21447
	},
	{
		city: 'Pratt',
		state: 'Kansas',
		lat: 37.646026,
		lon: -98.734756276566756,
		pop: 6835
	},
	{
		city: 'Prescott',
		state: 'Kansas',
		lat: 38.0631325,
		lon: -94.69539511337689,
		pop: 264
	},
	{
		city: 'Preston',
		state: 'Kansas',
		lat: 37.758264,
		lon: -98.555439460666662,
		pop: 158
	},
	{
		city: 'Pretty Prairie',
		state: 'Kansas',
		lat: 37.7786055,
		lon: -98.023463320973548,
		pop: 680
	},
	{
		city: 'Princeton',
		state: 'Kansas',
		lat: 38.4890195,
		lon: -95.271645287954556,
		pop: 277
	},
	{
		city: 'Protection',
		state: 'Kansas',
		lat: 37.201368,
		lon: -99.482056158845637,
		pop: 514
	},
	{
		city: 'Quenemo',
		state: 'Kansas',
		lat: 38.5793865,
		lon: -95.526453534626526,
		pop: 388
	},
	{
		city: 'Quinter',
		state: 'Kansas',
		lat: 39.063153,
		lon: -100.236753848330721,
		pop: 918
	},
	{
		city: 'Radium',
		state: 'Kansas',
		lat: 38.173683,
		lon: -98.893977397461924,
		pop: 25
	},
	{
		city: 'Ramona',
		state: 'Kansas',
		lat: 38.597924,
		lon: -97.063676678837552,
		pop: 187
	},
	{
		city: 'Randall',
		state: 'Kansas',
		lat: 39.6416545,
		lon: -98.045321731693917,
		pop: 65
	},
	{
		city: 'Randolph',
		state: 'Kansas',
		lat: 39.4289575,
		lon: -96.759728811529271,
		pop: 163
	},
	{
		city: 'Ransom',
		state: 'Kansas',
		lat: 38.635684,
		lon: -99.932572030868869,
		pop: 294
	},
	{
		city: 'Rantoul',
		state: 'Kansas',
		lat: 38.5489505,
		lon: -95.101052052140957,
		pop: 184
	},
	{
		city: 'Raymond',
		state: 'Kansas',
		lat: 38.278143,
		lon: -98.413566935752982,
		pop: 79
	},
	{
		city: 'Reading',
		state: 'Kansas',
		lat: 38.519295,
		lon: -95.95670254005168,
		pop: 231
	},
	{
		city: 'Redfield',
		state: 'Kansas',
		lat: 37.836233,
		lon: -94.881333853248265,
		pop: 146
	},
	{
		city: 'Republic',
		state: 'Kansas',
		lat: 39.924478,
		lon: -97.824943542418055,
		pop: 116
	},
	{
		city: 'Reserve',
		state: 'Kansas',
		lat: 39.976148,
		lon: -95.56505944468023,
		pop: 84
	},
	{
		city: 'Rexford',
		state: 'Kansas',
		lat: 39.4705105,
		lon: -100.743936188248512,
		pop: 232
	},
	{
		city: 'Richfield',
		state: 'Kansas',
		lat: 37.26529,
		lon: -101.782793050363836,
		pop: 43
	},
	{
		city: 'Richmond',
		state: 'Kansas',
		lat: 38.4012935,
		lon: -95.254923730807093,
		pop: 464
	},
	{
		city: 'Riley',
		state: 'Kansas',
		lat: 39.29881,
		lon: -96.828515848557686,
		pop: 939
	},
	{
		city: 'Riverton',
		state: 'Kansas',
		lat: 37.072473,
		lon: -94.70887975,
		pop: 929
	},
	{
		city: 'Robinson',
		state: 'Kansas',
		lat: 39.8156065,
		lon: -95.411786500727217,
		pop: 234
	},
	{
		city: 'Roeland Park',
		state: 'Kansas',
		lat: 39.033792,
		lon: -94.637762040269322,
		pop: 6731
	},
	{
		city: 'Rolla',
		state: 'Kansas',
		lat: 37.1182115,
		lon: -101.632566185957742,
		pop: 442
	},
	{
		city: 'Rosalia',
		state: 'Kansas',
		lat: 37.814466,
		lon: -96.625004844593548,
		pop: 171
	},
	{
		city: 'Rose Hill',
		state: 'Kansas',
		lat: 37.566878,
		lon: -97.136284968675767,
		pop: 3931
	},
	{
		city: 'Roseland',
		state: 'Kansas',
		lat: 37.2803,
		lon: -94.8441085,
		pop: 77
	},
	{
		city: 'Rossville',
		state: 'Kansas',
		lat: 39.136583,
		lon: -95.949737572704606,
		pop: 1151
	},
	{
		city: 'Roxbury',
		state: 'Kansas',
		lat: 38.550549,
		lon: -97.427412811931248,
		pop: 104
	},
	{
		city: 'Rozel',
		state: 'Kansas',
		lat: 38.19429,
		lon: -99.401925206497054,
		pop: 156
	},
	{
		city: 'Rush Center',
		state: 'Kansas',
		lat: 38.464227,
		lon: -99.309147959779551,
		pop: 170
	},
	{
		city: 'Russell',
		state: 'Kansas',
		lat: 38.8859175,
		lon: -98.853974740857936,
		pop: 24
	},
	{
		city: 'Russell Springs',
		state: 'Kansas',
		lat: 38.9127895,
		lon: -101.175855965765606,
		pop: 4506
	},
	{
		city: 'Sabetha',
		state: 'Kansas',
		lat: 39.900354,
		lon: -95.799706653279259,
		pop: 2571
	},
	{
		city: 'Salina',
		state: 'Kansas',
		lat: 38.8289435,
		lon: -97.603263996894412,
		pop: 47707
	},
	{
		city: 'Satanta',
		state: 'Kansas',
		lat: 37.4372855,
		lon: -100.987862647468702,
		pop: 1133
	},
	{
		city: 'Savonburg',
		state: 'Kansas',
		lat: 37.7495645,
		lon: -95.144408307083978,
		pop: 109
	},
	{
		city: 'Sawyer',
		state: 'Kansas',
		lat: 37.4971375,
		lon: -98.681361297026072,
		pop: 124
	},
	{
		city: 'Scammon',
		state: 'Kansas',
		lat: 37.278224,
		lon: -94.822085630174755,
		pop: 482
	},
	{
		city: 'Scandia',
		state: 'Kansas',
		lat: 39.7947225,
		lon: -97.782709767210136,
		pop: 372
	},
	{
		city: 'Schoenchen',
		state: 'Kansas',
		lat: 38.712646,
		lon: -99.331918157726847,
		pop: 207
	},
	{
		city: 'Scott City',
		state: 'Kansas',
		lat: 38.477964,
		lon: -100.901389846233755,
		pop: 3816
	},
	{
		city: 'Scottsville',
		state: 'Kansas',
		lat: 39.542568,
		lon: -97.952296721121101,
		pop: 25
	},
	{
		city: 'Scranton',
		state: 'Kansas',
		lat: 38.7771155,
		lon: -95.740024185767098,
		pop: 710
	},
	{
		city: 'Sedan',
		state: 'Kansas',
		lat: 37.128801,
		lon: -96.183283412160478,
		pop: 1124
	},
	{
		city: 'Sedgwick',
		state: 'Kansas',
		lat: 37.915171,
		lon: -97.420281020925643,
		pop: 1695
	},
	{
		city: 'Selden',
		state: 'Kansas',
		lat: 39.5416185,
		lon: -100.566833737284611,
		pop: 219
	},
	{
		city: 'Seneca',
		state: 'Kansas',
		lat: 39.836688,
		lon: -96.066233676630347,
		pop: 1991
	},
	{
		city: 'Severance',
		state: 'Kansas',
		lat: 39.7664505,
		lon: -95.249312993018606,
		pop: 94
	},
	{
		city: 'Severy',
		state: 'Kansas',
		lat: 37.621804,
		lon: -96.229437848226951,
		pop: 259
	},
	{
		city: 'Seward',
		state: 'Kansas',
		lat: 38.179312,
		lon: -98.794354913910212,
		pop: 64
	},
	{
		city: 'Sharon',
		state: 'Kansas',
		lat: 37.250327,
		lon: -98.418364230547539,
		pop: 748
	},
	{
		city: 'Sharon Springs',
		state: 'Kansas',
		lat: 38.895219,
		lon: -101.752485995448637,
		pop: 158
	},
	{
		city: 'Shawnee',
		state: 'Kansas',
		lat: 39.015962,
		lon: -94.796877144938094,
		pop: 62209
	},
	{
		city: 'Silver Lake',
		state: 'Kansas',
		lat: 39.099014,
		lon: -95.854511194053558,
		pop: 1439
	},
	{
		city: 'Simpson',
		state: 'Kansas',
		lat: 39.386124,
		lon: -97.934168802677021,
		pop: 86
	},
	{
		city: 'Smith Center',
		state: 'Kansas',
		lat: 39.779748,
		lon: -98.780994169829114,
		pop: 1665
	},
	{
		city: 'Smolan',
		state: 'Kansas',
		lat: 38.737706,
		lon: -97.684148307458969,
		pop: 215
	},
	{
		city: 'Soldier',
		state: 'Kansas',
		lat: 39.536879,
		lon: -95.964973021629689,
		pop: 136
	},
	{
		city: 'Solomon',
		state: 'Kansas',
		lat: 38.9211155,
		lon: -97.372097235371285,
		pop: 1095
	},
	{
		city: 'South Haven',
		state: 'Kansas',
		lat: 37.0527555,
		lon: -97.400943706458463,
		pop: 363
	},
	{
		city: 'South Hutchinson',
		state: 'Kansas',
		lat: 38.0242225,
		lon: -97.938191566224248,
		pop: 2457
	},
	{
		city: 'Spearville',
		state: 'Kansas',
		lat: 37.8476275,
		lon: -99.755360928539176,
		pop: 773
	},
	{
		city: 'Speed',
		state: 'Kansas',
		lat: 39.67684,
		lon: -99.420042193068554,
		pop: 37
	},
	{
		city: 'Spivey',
		state: 'Kansas',
		lat: 37.448415,
		lon: -98.165716714628857,
		pop: 78
	},
	{
		city: 'Spring Hill',
		state: 'Kansas',
		lat: 38.756421,
		lon: -94.836925530339101,
		pop: 5437
	},
	{
		city: 'St. Francis',
		state: 'Kansas',
		lat: 39.771062,
		lon: -101.802098969827597,
		pop: 1329
	},
	{
		city: 'St. George',
		state: 'Kansas',
		lat: 39.191732,
		lon: -96.416082040892192,
		pop: 639
	},
	{
		city: 'St. John',
		state: 'Kansas',
		lat: 37.999511,
		lon: -98.761320401830588,
		pop: 1295
	},
	{
		city: 'St. Marys',
		state: 'Kansas',
		lat: 39.195727,
		lon: -96.071731617152153,
		pop: 2627
	},
	{
		city: 'St. Paul',
		state: 'Kansas',
		lat: 37.5195075,
		lon: -95.17605009715777,
		pop: 629
	},
	{
		city: 'Stafford',
		state: 'Kansas',
		lat: 37.961221,
		lon: -98.599578431969263,
		pop: 1042
	},
	{
		city: 'Stark',
		state: 'Kansas',
		lat: 37.6891275,
		lon: -95.1428878736846,
		pop: 72
	},
	{
		city: 'Sterling',
		state: 'Kansas',
		lat: 38.2103595,
		lon: -98.205512679863489,
		pop: 2328
	},
	{
		city: 'Stockton',
		state: 'Kansas',
		lat: 39.436512,
		lon: -99.270694765184373,
		pop: 1329
	},
	{
		city: 'Strong City',
		state: 'Kansas',
		lat: 38.3949225,
		lon: -96.535953394736822,
		pop: 485
	},
	{
		city: 'Sublette',
		state: 'Kansas',
		lat: 37.4822145,
		lon: -100.844200276395981,
		pop: 1453
	},
	{
		city: 'Summerfield',
		state: 'Kansas',
		lat: 39.9949415,
		lon: -96.350679944043122,
		pop: 156
	},
	{
		city: 'Sun City',
		state: 'Kansas',
		lat: 37.378522,
		lon: -98.916448022624422,
		pop: 53
	},
	{
		city: 'Susank',
		state: 'Kansas',
		lat: 38.6403155,
		lon: -98.774529446880962,
		pop: 34
	},
	{
		city: 'Sylvan Grove',
		state: 'Kansas',
		lat: 39.011451,
		lon: -98.392581788708299,
		pop: 279
	},
	{
		city: 'Sylvia',
		state: 'Kansas',
		lat: 37.958738,
		lon: -98.409251074786312,
		pop: 218
	},
	{
		city: 'Syracuse',
		state: 'Kansas',
		lat: 37.9551705,
		lon: -101.79806299365967,
		pop: 1812
	},
	{
		city: 'Talmage',
		state: 'Kansas',
		lat: 39.0267755,
		lon: -97.259647381578944,
		pop: 99
	},
	{
		city: 'Tampa',
		state: 'Kansas',
		lat: 38.547208,
		lon: -97.154311145965863,
		pop: 112
	},
	{
		city: 'Tescott',
		state: 'Kansas',
		lat: 39.0144835,
		lon: -97.877856624738669,
		pop: 319
	},
	{
		city: 'Thayer',
		state: 'Kansas',
		lat: 37.481346,
		lon: -95.495800210175517,
		pop: 497
	},
	{
		city: 'Timken',
		state: 'Kansas',
		lat: 38.472104,
		lon: -99.178983795744259,
		pop: 76
	},
	{
		city: 'Tipton',
		state: 'Kansas',
		lat: 39.339429,
		lon: -98.470829273429615,
		pop: 210
	},
	{
		city: 'Tonganoxie',
		state: 'Kansas',
		lat: 39.1116805,
		lon: -95.0826342156995,
		pop: 4996
	},
	{
		city: 'Topeka',
		state: 'Kansas',
		lat: 39.036832,
		lon: -95.692601939328981,
		pop: 127473
	},
	{
		city: 'Toronto',
		state: 'Kansas',
		lat: 37.7982475,
		lon: -95.949994244212093,
		pop: 281
	},
	{
		city: 'Towanda',
		state: 'Kansas',
		lat: 37.792278,
		lon: -97.001808006689089,
		pop: 1450
	},
	{
		city: 'Treece',
		state: 'Kansas',
		lat: 37.0002965,
		lon: -94.84373975,
		pop: 138
	},
	{
		city: 'Tribune',
		state: 'Kansas',
		lat: 38.4723055,
		lon: -101.754153190138112,
		pop: 741
	},
	{
		city: 'Troy',
		state: 'Kansas',
		lat: 39.781973,
		lon: -95.102621095122629,
		pop: 1010
	},
	{
		city: 'Turon',
		state: 'Kansas',
		lat: 37.807849,
		lon: -98.42864447402772,
		pop: 387
	},
	{
		city: 'Tyro',
		state: 'Kansas',
		lat: 37.037911,
		lon: -95.82174676595244,
		pop: 220
	},
	{
		city: 'Udall',
		state: 'Kansas',
		lat: 37.3910405,
		lon: -97.114677235721985,
		pop: 746
	},
	{
		city: 'Ulysses',
		state: 'Kansas',
		lat: 37.5740265,
		lon: -101.357807873185379,
		pop: 6161
	},
	{
		city: 'Uniontown',
		state: 'Kansas',
		lat: 37.8472955,
		lon: -94.974509778427105,
		pop: 272
	},
	{
		city: 'Utica',
		state: 'Kansas',
		lat: 38.6426785,
		lon: -100.170002408586441,
		pop: 158
	},
	{
		city: 'Valley Center',
		state: 'Kansas',
		lat: 37.828683,
		lon: -97.379947851551478,
		pop: 6822
	},
	{
		city: 'Valley Falls',
		state: 'Kansas',
		lat: 39.3399415,
		lon: -95.461971743009656,
		pop: 1192
	},
	{
		city: 'Vassar',
		state: 'Kansas',
		lat: 38.6512085,
		lon: -95.607393426056348,
		pop: 530
	},
	{
		city: 'Vermillion',
		state: 'Kansas',
		lat: 39.718309,
		lon: -96.266450012674824,
		pop: 112
	},
	{
		city: 'Victoria',
		state: 'Kansas',
		lat: 38.8526725,
		lon: -99.146366400570486,
		pop: 1214
	},
	{
		city: 'Vining',
		state: 'Kansas',
		lat: 39.5680785,
		lon: -97.293632363836593,
		pop: 45
	},
	{
		city: 'Viola',
		state: 'Kansas',
		lat: 37.4827605,
		lon: -97.644527226800946,
		pop: 130
	},
	{
		city: 'Virgil',
		state: 'Kansas',
		lat: 37.9803535,
		lon: -96.011619571816638,
		pop: 71
	},
	{
		city: 'WaKeeney',
		state: 'Kansas',
		lat: 39.0227885,
		lon: -99.879365579264586,
		pop: 1862
	},
	{
		city: 'Wakarusa',
		state: 'Kansas',
		lat: 38.891647,
		lon: -95.699433381082059,
		pop: 260
	},
	{
		city: 'Wakefield',
		state: 'Kansas',
		lat: 39.2160435,
		lon: -97.016091856339401,
		pop: 980
	},
	{
		city: 'Waldo',
		state: 'Kansas',
		lat: 39.120024,
		lon: -98.797970344980939,
		pop: 30
	},
	{
		city: 'Waldron',
		state: 'Kansas',
		lat: 37.002692,
		lon: -98.182032570700613,
		pop: 11
	},
	{
		city: 'Wallace',
		state: 'Kansas',
		lat: 38.912475,
		lon: -101.593014857665636,
		pop: 57
	},
	{
		city: 'Walnut',
		state: 'Kansas',
		lat: 37.599614,
		lon: -95.074741743989605,
		pop: 220
	},
	{
		city: 'Walton',
		state: 'Kansas',
		lat: 38.119996,
		lon: -97.256308284161122,
		pop: 235
	},
	{
		city: 'Wamego',
		state: 'Kansas',
		lat: 39.205122,
		lon: -96.310511607287452,
		pop: 4372
	},
	{
		city: 'Washington',
		state: 'Kansas',
		lat: 39.8175845,
		lon: -97.054006123987563,
		pop: 1131
	},
	{
		city: 'Waterville',
		state: 'Kansas',
		lat: 39.690956,
		lon: -96.748570349202822,
		pop: 680
	},
	{
		city: 'Wathena',
		state: 'Kansas',
		lat: 39.7613035,
		lon: -94.94259330378236,
		pop: 1364
	},
	{
		city: 'Waverly',
		state: 'Kansas',
		lat: 38.395795,
		lon: -95.605358151324168,
		pop: 592
	},
	{
		city: 'Webber',
		state: 'Kansas',
		lat: 39.9342495,
		lon: -98.035276782611533,
		pop: 25
	},
	{
		city: 'Weir',
		state: 'Kansas',
		lat: 37.3071135,
		lon: -94.773490957972342,
		pop: 686
	},
	{
		city: 'Welda',
		state: 'Kansas',
		lat: 38.174311,
		lon: -95.289550883001041,
		pop: 129
	},
	{
		city: 'Wellington',
		state: 'Kansas',
		lat: 37.274727,
		lon: -97.398314603263657,
		pop: 8172
	},
	{
		city: 'Wellsville',
		state: 'Kansas',
		lat: 38.716202,
		lon: -95.07984556957102,
		pop: 1857
	},
	{
		city: 'Weskan',
		state: 'Kansas',
		lat: 38.8660775,
		lon: -101.968486248057616,
		pop: 161
	},
	{
		city: 'West Mineral',
		state: 'Kansas',
		lat: 37.285581,
		lon: -94.930156094450226,
		pop: 185
	},
	{
		city: 'Westmoreland',
		state: 'Kansas',
		lat: 39.3930355,
		lon: -96.412747166588844,
		pop: 778
	},
	{
		city: 'Westphalia',
		state: 'Kansas',
		lat: 38.182514,
		lon: -95.49124385869564,
		pop: 163
	},
	{
		city: 'Westwood',
		state: 'Kansas',
		lat: 39.038788,
		lon: -94.616557935651684,
		pop: 359
	},
	{
		city: 'Westwood Hills',
		state: 'Kansas',
		lat: 39.038788,
		lon: -94.609662813118803,
		pop: 1506
	},
	{
		city: 'Wetmore',
		state: 'Kansas',
		lat: 39.635723,
		lon: -95.811978440751446,
		pop: 368
	},
	{
		city: 'Wheaton',
		state: 'Kansas',
		lat: 39.502162,
		lon: -96.318985375635805,
		pop: 95
	},
	{
		city: 'White City',
		state: 'Kansas',
		lat: 38.792757,
		lon: -96.73500533922261,
		pop: 618
	},
	{
		city: 'White Cloud',
		state: 'Kansas',
		lat: 39.975876,
		lon: -95.300116155172418,
		pop: 176
	},
	{
		city: 'Whitewater',
		state: 'Kansas',
		lat: 37.9624565,
		lon: -97.145764884152101,
		pop: 718
	},
	{
		city: 'Whiting',
		state: 'Kansas',
		lat: 39.588636,
		lon: -95.611147440615781,
		pop: 187
	},
	{
		city: 'Wichita',
		state: 'Kansas',
		lat: 37.681309,
		lon: -97.387208902503801,
		pop: 382368
	},
	{
		city: 'Willard',
		state: 'Kansas',
		lat: 39.0943135,
		lon: -95.942764884458668,
		pop: 92
	},
	{
		city: 'Williamsburg',
		state: 'Kansas',
		lat: 38.4819525,
		lon: -95.470780047140579,
		pop: 397
	},
	{
		city: 'Willis',
		state: 'Kansas',
		lat: 39.722524,
		lon: -95.505936666666656,
		pop: 38
	},
	{
		city: 'Willowbrook',
		state: 'Kansas',
		lat: 38.1011995,
		lon: -97.991808384433952,
		pop: 87
	},
	{
		city: 'Wilmore',
		state: 'Kansas',
		lat: 37.3346135,
		lon: -99.209110278729568,
		pop: 53
	},
	{
		city: 'Wilroads Gardens',
		state: 'Kansas',
		lat: 37.715109,
		lon: -99.923442284976673,
		pop: 609
	},
	{
		city: 'Wilsey',
		state: 'Kansas',
		lat: 38.6356075,
		lon: -96.675922165238092,
		pop: 153
	},
	{
		city: 'Wilson',
		state: 'Kansas',
		lat: 38.826128,
		lon: -98.474682891169152,
		pop: 781
	},
	{
		city: 'Winchester',
		state: 'Kansas',
		lat: 39.3229565,
		lon: -95.268815102459058,
		pop: 551
	},
	{
		city: 'Windom',
		state: 'Kansas',
		lat: 38.3837525,
		lon: -97.910942906224278,
		pop: 130
	},
	{
		city: 'Winfield',
		state: 'Kansas',
		lat: 37.2354365,
		lon: -96.986961992774567,
		pop: 12301
	},
	{
		city: 'Winona',
		state: 'Kansas',
		lat: 39.0619495,
		lon: -101.243801506723713,
		pop: 162
	},
	{
		city: 'Woodbine',
		state: 'Kansas',
		lat: 38.7953645,
		lon: -96.95958525,
		pop: 170
	},
	{
		city: 'Woodston',
		state: 'Kansas',
		lat: 39.4534705,
		lon: -99.098307160333704,
		pop: 136
	},
	{
		city: 'Wright',
		state: 'Kansas',
		lat: 37.774884,
		lon: -99.890806173188651,
		pop: 163
	},
	{
		city: 'Yates Center',
		state: 'Kansas',
		lat: 37.878447,
		lon: -95.732115106094398,
		pop: 1417
	},
	{
		city: 'Yoder',
		state: 'Kansas',
		lat: 37.9440445,
		lon: -97.866367645284328,
		pop: 194
	},
	{
		city: 'Zenda',
		state: 'Kansas',
		lat: 37.4442695,
		lon: -98.280510493987407,
		pop: 90
	},
	{
		city: 'Zurich',
		state: 'Kansas',
		lat: 39.231236,
		lon: -99.434717149559503,
		pop: 99
	}
];

var field = VectorField.read(windData, true);

var mapAnimator;
var magicNumber = 3.5;
magicNumber = 2.23693629;
var legendSpeeds = [5 / magicNumber, 10 / magicNumber, 15 / magicNumber, 20 / magicNumber, 25 / magicNumber, 30 / magicNumber];

var MapMask = function(image, width, height) {
	this.image = image;
	this.width = width;
	this.height = height;
};

MapMask.prototype.endMove = function(animator) {
	this.move(animator);
};

MapMask.prototype.move = function(animator) {
	var s = this.image.style;
	s.width = ~~(animator.scale * this.width) + 'px';
	s.height = ~~(animator.scale * this.height) + 'px';
	s.left = animator.dx + 'px';
	s.top = animator.dy + 'px';
};

function isAnimating() {
	return document.getElementById('animating').checked;
}

function showCities() {
	document.getElementById('city-display').style.visibility = document.getElementById('show-cities').checked ? 'visible' : 'hidden';
}

function doUnzoom() {
	mapAnimator.unzoom();
}

function format(x) {
	x = Math.round(x * 10) / 10;
	var a1 = ~~x;
	var a2 = (~~(x * 10)) % 10;
	return a1 + '.' + a2;
}

function init() {
	loading = false;
	var timestamp = windData.timestamp || 'unknown on unknown';
	var parts = timestamp.split('on');
	var time = parts[0].trim();
	var day = parts[1].trim().replace(' 0', ' ');
	day = day.replace('September', 'Sept.');
	day = day.replace('November', 'Nov.');
	day = day.replace('December', 'Dec.');
	
	document.getElementById('update-time').innerHTML = '<span id="day">' + day + '</span><br>' + time + ' CST' + '<br><span id="time-explanation">(time of forecast download)</span>';
	
	var avg = field.averageLength * 2.23693629;
	var max = field.maxLength * 2.23693629;
	
	document.getElementById('average-speed').innerHTML = '<br>top speed: <b>' + format(max) + ' mph</b><br>' + 'average: <b>' + format(avg) + ' mph</b>';

	var canvas = document.getElementById('display');
	var imageCanvas = document.getElementById('image-canvas');
	var mapProjection = new ScaledAlbers(8850, 0, canvas.height - 8, -102, 36.99);

	var isMacFF = navigator.platform.indexOf('Mac') != -1 && navigator.userAgent.indexOf('Firefox') != -1;
	var isWinFF = navigator.platform.indexOf('Win') != -1 && navigator.userAgent.indexOf('Firefox') != -1;
	var isWinIE = navigator.platform.indexOf('Win') != -1 && navigator.userAgent.indexOf('MSIE') != -1;

	var numParticles = isMacFF || isWinIE ? 3500 : 5000;
	var display = new MotionDisplay(canvas, imageCanvas, field, numParticles, mapProjection);

	if (isWinFF || isWinIE) {
		display.setAlpha(0.05);
	}

	var navDiv = document.getElementById('city-display');
	var unzoom = document.getElementById('unzoom');
	mapAnimator = new Animator(navDiv, isAnimating, unzoom);
	mapAnimator.add(display);

	var mask = new MapMask(document.getElementById('mask'), 900, 461);
	mapAnimator.add(mask);

	var callout = document.getElementById('callout');
	var hovercard = new MotionDetails(navDiv, callout, field, mapProjection, mapAnimator);

	var cityCanvas = document.getElementById('city-display');
	cityDisplay = new CityDisplay(cities, cityCanvas, mapProjection);
	mapAnimator.add(cityDisplay);
	cityDisplay.move();

	var legendAnimator = new Animator(null, isAnimating);

	var speedScaleFactor = 20 * 2.23693629;
	speedScaleFactor = 10 * 2.23693629;
	for (var i = 1; i <= legendSpeeds.length; i++) {
		var c = document.getElementById('legend' + i);
		var legendField = VectorField.constant(legendSpeeds[i - 1] * speedScaleFactor, 0, 0, 0, c.width, c.height);
		var legend = new MotionDisplay(c, null, legendField, 30);
		legend.maxLength = field.maxLength * speedScaleFactor;
		legendAnimator.add(legend);
	}

	mapAnimator.start(20);
	legendAnimator.start(20);
}