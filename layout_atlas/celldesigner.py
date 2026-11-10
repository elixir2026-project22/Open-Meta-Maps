"""Conservative CellDesigner center-line decoder; immutable observations retained.

Non-branch direct reactions only. Geometry convention checked against the
CellDesigner4 extension specification and cd2sbgnml GeometryUtils. Unsupported
branch/anchor clipping is explicitly excluded, not replaced with invented paths.
"""
import xml.etree.ElementTree as ET
from collections import Counter
from core import *
import master_map_audit


def local_to_global(start, end, points):
    a=np.array(start,float);v=np.array(end,float)-a
    if np.linalg.norm(v)<1e-8:raise ValueError('degenerate_endpoint_frame')
    perpendicular=np.array([-v[1],v[0]])
    return [(a+float(x)*v+float(y)*perpendicular).tolist() for x,y in points]


def decode(path):
    raw=master_map_audit.parse(path);root=ET.parse(path).getroot()
    records={r['id']:r for r in raw['reactions']};audit=Counter();events=[];glyphs={}
    for aid,alias in sorted(raw['aliases'].items()):
        s=raw['species'][alias['species']];b=alias['bounds']
        glyphs[aid]=dict(key=aid,entity=s['identity'],label=s['name'],
                         xy=[b['x']+b['w']/2,b['y']+b['h']/2],bounds=b,
                         source_species=alias['species'],observed=True,position_origin='observed')
    excluded=[]
    for e in root.iter():
        if local(e)!='reaction' or not e.tag.startswith('{http://www.sbml.org/sbml/'):continue
        r=records[e.get('id')];ext=next((x for x in e.iter() if local(x)=='extension'),None)
        if ext is None:excluded.append(dict(id=r['id'],reason='missing_extension'));continue
        children={local(x):x for x in ext}
        reactants=[x for x in ext.iter() if local(x)=='baseReactant']
        products=[x for x in ext.iter() if local(x)=='baseProduct']
        reason=None
        if (len(reactants),len(products))!=(1,1):reason='branch_geometry_not_decoded'
        elif not r['equation_numeric_complete']:reason='incomplete_participant_identity_or_stoichiometry'
        elif not all(x['valid'] for x in r['glyph_links']):reason='invalid_alias_link'
        else:
            scheme=children.get('connectScheme')
            if scheme is not None and scheme.get('connectPolicy')!='direct':reason='non_direct_scheme'
        if reason:excluded.append(dict(id=r['id'],reason=reason));audit[reason]+=1;continue
        a=reactants[0].get('alias');b=products[0].get('alias')
        try:
            ep=children.get('editPoints');text=(ep.text or '') if ep is not None else ''
            if ep is not None and any(k in ep.attrib for k in ('num0','num1','num2','tShapeIndex')):raise ValueError('branch_attributes_on_nonbranch')
            points=[list(map(float,p.split(','))) for p in text.split()]
            if any(len(p)!=2 for p in points):raise ValueError('invalid_control_point')
            interior=local_to_global(glyphs[a]['xy'],glyphs[b]['xy'],points)
            path_points=[glyphs[a]['xy'],*interior,glyphs[b]['xy']]
            index=int(scheme.get('rectangleIndex','0')) if scheme is not None else 0
            if not 0<=index<len(path_points)-1:raise ValueError('invalid_process_segment')
            center=((np.array(path_points[index])+path_points[index+1])/2).tolist()
            if not np.isfinite(path_points).all():raise ValueError('nonfinite_route')
        except (ValueError,KeyError) as err:
            reason=str(err);excluded.append(dict(id=r['id'],reason=reason));audit[reason]+=1;continue
        visible={}
        for link in r['glyph_links']:
            if link['kind'] not in ('baseReactant','baseProduct'):continue
            entity=glyphs[link['alias']]['entity'];coefficient=r['coefficients'].get(entity)
            if coefficient is None:reason='core_participant_cancels_or_unresolved';break
            visible[link['alias']]=coefficient
        if reason:excluded.append(dict(id=r['id'],reason=reason));audit[reason]+=1;continue
        events.append(dict(key='r:'+r['id'],id=r['name'],source_reaction_id=r['id'],
             coefficients=r['coefficients'],references={'vmh.label':[r['name']]},
             reversible=r['reversible'],reversibility_origin=r['reversibility_origin'],
             visible=visible,xy=center,position_origin='decoded_process_segment_midpoint_not_observed',
             core=list(visible),path_centerline=path_points,process_segment=index,
             link_anchor_observations=[dict(kind=local(x),alias=x.get('alias'),anchors=[dict(ch.attrib) for ch in x if local(ch)=='linkAnchor']) for x in reactants+products],
             all_glyph_links=r['glyph_links'],raw_edit_points=r['paths'],
             route_scope='source center-line controls; glyph-boundary clipping and side-route decoding not implemented'))
        audit['decoded_with_controls' if points else 'decoded_direct_no_controls']+=1
    return dict(id=uid(['micromap',sha(Path(path).read_bytes())]),source='micromap',
       species='microbial master design; not one organism',name='MicroMap master v4.0',
       lineage='micromap:'+sha(Path(path).read_bytes()),source_sha256=sha(Path(path).read_bytes()),
       path=str(Path(path).relative_to(BASE)),license='CC-BY-NC-SA-4.0',
       partition='research_extension',admitted_to_open_core=False,
       glyphs=glyphs,events=events,raw_reactions=raw['reactions'],
       audit=dict(audit),excluded=excluded,raw_species=raw['species'])
