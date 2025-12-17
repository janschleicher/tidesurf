{{ fullname | escape | underline}}
{% set non_init_methods = methods | reject("equalto", "__init__") | list %}
.. currentmodule:: {{ module }}

.. add toctree option to make autodoc generate the pages

.. autoclass:: {{ objname }}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Attributes

   .. autosummary::
      :toctree: .
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block methods %}
   {% if non_init_methods %}
   .. rubric:: Methods

   .. autosummary::
      :toctree: .
   {% for item in methods %}
      {%- if item != '__init__' %}
      ~{{ name }}.{{ item }}
      {%- endif -%}
   {%- endfor %}
   {% endif %}
   {% endblock %}