{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
    :inherited-members:
    :show-inheritance:
    :members:

    {% if methods %}
    .. rubric:: Methods
    
    .. autosummary::
    {% for item in methods %}
        ~{{ name}}.{{ item }}
    {%- endfor %}
    {% endif %}
    .. automethod:: __init__ 
