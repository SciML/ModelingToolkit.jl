block SimpleAdder
	parameter Real value = 1.0;
	input Real a;
	input Real b;
	Real c(start = 1.0, fixed = true);
	output Real out;
	output Real out2;
equation
	der(c) = out;
	out = a + b + value;
	out2 = 2 * c;
end SimpleAdder;
