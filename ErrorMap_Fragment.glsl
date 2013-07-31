#extension GL_EXT_gpu_shader4 : enable

varying vec3 fragment_normal;
//uniform vec3 L0_pos;
uniform vec3 L1_pos;
uniform sampler1D colorMap;

vec4 HSLtoRGB(vec4 HSL)
{
	float H = HSL.x;
	float S = HSL.y;
	float L = HSL.z;

	float temp2;
	if(L < 0.5)
		temp2 = L*(1.0+S);
	else
		temp2 = L+S - L*S;

	float temp1 = 2.0*L - temp2;

	vec3 temp3 = vec3(H+1.0/3.0, H, H-1.0/3.0);
	if(temp3.r < 0.0) temp3.r = temp3.r+1.0;
	if(temp3.g < 0.0) temp3.g = temp3.g+1.0;
	if(temp3.b < 0.0) temp3.b = temp3.b+1.0;

	if(temp3.r > 1.0) temp3.r = temp3.r - 1.0;
	if(temp3.g > 1.0) temp3.g = temp3.g - 1.0;
	if(temp3.b > 1.0) temp3.b = temp3.b - 1.0;

	vec4 result;
	if(6.0*temp3.r < 1.0) result.r = temp1 +(temp2 - temp1)*6.0*temp3.r;
	else if(2.0*temp3.r < 1.0) result.r = temp2;
	else if(3.0*temp3.r < 2.0) result.r = temp1+(temp2-temp1)*((2.0/3.0) - temp3.r)*6.0;
	else result.r = temp1;

	if(6.0*temp3.g < 1.0) result.g = temp1 +(temp2 - temp1)*6.0*temp3.g;
	else if(2.0*temp3.g < 1.0) result.g = temp2;
	else if(3.0*temp3.g < 2.0) result.g = temp1+(temp2-temp1)*((2.0/3.0) - temp3.g)*6.0;
	else result.g = temp1;

	if(6.0*temp3.b < 1.0) result.b = temp1 +(temp2 - temp1)*6.0*temp3.b;
	else if(2.0*temp3.b < 1.0) result.b = temp2;
	else if(3.0*temp3.b < 2.0) result.b = temp1+(temp2-temp1)*((2.0/3.0) - temp3.b)*6.0;
	else result.b = temp1;

	result.a = 0.0;
	return result;
}

float LightIntensity()
{
	//vec3 L0_pos = vec3(1.0, 0.0, 0.0);
	//vec3 L1_pos = vec3(0.0, 0.0, 1.0);
	//float L0 = max(dot(fragment_normal, L0_pos), 0.0);
	float L1 = max(dot(fragment_normal, L1_pos), 0.0);

	//float total = L0 + L1;
	//if(total > 1.0)
	//	total = 1.0;

	return L1;
	//return 0.5;

}

void main(void)
{
	float error = gl_TexCoord[0].x;

	float light = LightIntensity();

	gl_FragColor = texture1D(colorMap, error)*light;
}
